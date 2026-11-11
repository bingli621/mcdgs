"""The neutrons that reach the sample of T-Rex, from ESS_butterfly and from Polygon_ESS_butterfly.

Both runs are :class:`.TRexToSample`, T-Rex up to its sample monitor, as in
:mod:`.diagnostics.trex.neutrons_on_sample`, with the same number of rays and the same seed. The first
has the source of :func:`.sources.source`, ``ESS_butterfly``. The second has that of
:func:`.sources.masked_source`, ``Polygon_ESS_butterfly``, which emits only where the six
choppers let neutrons through and multiplies each ray's weight by the acceptance, the share of
the sampled (inverse velocity, emission time) plane the choppers leave open. Both should give the
sample the same intensity, the masked source from many more rays that reach it.

The sample monitor of both has wavelength bins of ``LAMBDA_RESOLUTION`` and time bins of
``TIME_RESOLUTION``, finer than those of T-Rex, to show the shape of each sub-pulse.
:func:`find_peaks` finds the sub-pulses in the time profile, and :func:`plot_sources` draws the
wavelength spectrum and the time profile of both, on a linear and on a log intensity axis,
``LOG_DECADES`` deep, and :func:`plot_peaks` the profile around each sub-pulse,
``PEAK_HALF_WIDTH`` on each side. :func:`summary_table` lists, for each source, the total at
the sample monitor and its relative error, the time of the run and the figure of merit, one over
the relative error squared times the time: how much a second of each run is worth.

:func:`main` runs both and saves the figures and the table. Run the module to do it::

    pixi run python -m mcdgs.diagnostics.trex.masked_source
"""

import json
import math
import time
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import mcstasscript as ms
import numpy as np
import scipp as sc
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from numpy.typing import NDArray

from mcdgs.diagnostics.cli import markdown_table, output_folder, parse_out, save_figure, save_text
from mcdgs.diagnostics.compare import read_mcstas_map
from mcdgs.diagnostics.trex.neutrons_on_sample import TRexToSample
from mcdgs.instruments.trex.monitors import sample_monitor
from mcdgs.instruments.trex.sources import wavelength_band
from mcdgs.mcstas.simulation import Simulation
from mcdgs.tof.reader import monitor_bins

MONITOR = "sample_monitor"  # where both sources are compared
LAMBDA_RESOLUTION = 0.002  # Å, the wavelength bin of its map: a tenth of that of T-Rex
TIME_RESOLUTION = 1e-5  # s, its time bin: a fifth of that of T-Rex
PEAK_THRESHOLD = 0.01  # a sub-pulse has bins above this fraction of the highest bin of the profile
PEAK_GAP = 1.0  # ms, a gap longer than this between such bins starts the next sub-pulse
PEAK_HALF_WIDTH = 0.1  # ms, the zoom on a sub-pulse shows this much on each side of its peak
PEAK_COLUMNS = 6  # the zooms on the sub-pulses are drawn this many to a row
LOG_DECADES = 6  # a log axis spans this many decades below the highest point


@dataclass(frozen=True)
class SourceRun:
    """The run of one source, at the sample monitor.

    Parameters
    ----------
    component:
        The McStas component of the source.
    folder:
        The run folder.
    seconds:
        The wall time of the run, compiling included, in s.
    intensity:
        The intensity per bin of the sample monitor, in n/s, of shape (time bins, wavelength
        bins).
    error:
        Its error per bin.
    acceptance:
        The acceptance of the masked source, from its region file, or ``None`` for the plain one.
    """

    component: str
    folder: Path
    seconds: float
    intensity: NDArray
    error: NDArray
    acceptance: float | None

    @property
    def total(self) -> float:
        """The total intensity at the sample monitor, in n/s."""
        return float(self.intensity.sum())

    @property
    def relative_error(self) -> float:
        """The relative error of :attr:`total`: the errors of the bins add in quadrature."""
        return float(np.sqrt((self.error**2).sum()) / self.total)

    @property
    def figure_of_merit(self) -> float:
        """One over the relative error squared times the time, in 1/s."""
        return 1 / (self.relative_error**2 * self.seconds)


def run_source(instrument: TRexToSample, folder_name: str, ncount: float, seed: int) -> SourceRun:
    """Run one source to the sample monitor, and read the monitor.

    Parameters
    ----------
    instrument:
        T-Rex to the sample monitor, with the source to run.
    folder_name:
        The run folder name.
    ncount:
        Number of McStas rays.
    seed:
        Seed of the random numbers.

    Returns
    -------
    :
        The run.
    """
    start = time.perf_counter()
    folder = Simulation(instrument).run(folder_name, ncount=ncount, seed=seed).parent
    seconds = time.perf_counter() - start
    intensity, error = read_mcstas_map(
        MONITOR, ms.load_data(str(folder)), monitor_bins(instrument.sample_monitor)
    )
    region = folder / f"{instrument.source.name}.json"  # the masked source writes its region
    acceptance = json.loads(region.read_text())["acceptance"] if region.is_file() else None
    return SourceRun(
        instrument.source.component_name, folder, seconds, intensity, error, acceptance
    )


def summary_table(plain: SourceRun, masked: SourceRun, ncount: float) -> str:
    """Return the totals, errors, times and figures of merit of both runs as a Markdown table.

    Parameters
    ----------
    plain, masked:
        The runs of both sources.
    ncount:
        Number of McStas rays of each.

    Returns
    -------
    :
        The table, a row per source, then the ratio of the masked to the plain one.
    """
    header = ["Source", "Rays", "At the sample", "Relative error", "Time", "Figure of merit"]
    rows = [
        (
            f"`{run.component}`"
            + ("" if run.acceptance is None else f", acceptance {run.acceptance:.4g}"),
            f"{ncount:.0e}",
            f"{run.total:.4g} n/s",
            f"{run.relative_error:.2%}",
            f"{run.seconds:.0f} s",
            f"{run.figure_of_merit:.3g} /s",
        )
        for run in (plain, masked)
    ]
    rows.append(
        (
            "Masked over plain",
            "",
            f"{masked.total / plain.total:.4f}",
            f"{masked.relative_error / plain.relative_error:.3g}",
            f"{masked.seconds / plain.seconds:.3g}",
            f"{masked.figure_of_merit / plain.figure_of_merit:.3g}",
        )
    )
    return markdown_table(header, rows)


def _profiles(
    bins: Mapping[str, sc.Variable], runs: tuple[SourceRun, ...]
) -> tuple[NDArray, list[tuple[NDArray, NDArray]]]:
    """Return the centers of the time bins, in ms, and the time profile of each run, with its error.

    Parameters
    ----------
    bins:
        The bin edges of the maps, ``"toa"`` and ``"wavelength"``.
    runs:
        The runs.

    Returns
    -------
    :
        The centers, and for each run its intensity per time bin, summed over wavelength, and
        its error: the errors of the bins add in quadrature.
    """
    times = sc.midpoints(bins["toa"]).to(unit="ms").values
    return times, [(run.intensity.sum(axis=1), np.sqrt((run.error**2).sum(axis=1))) for run in runs]


def find_peaks(
    run: SourceRun,
    bins: Mapping[str, sc.Variable],
    threshold: float = PEAK_THRESHOLD,
    gap: float = PEAK_GAP,
) -> list[float]:
    """Find the sub-pulses in the time profile of a run.

    A sub-pulse is a run of time bins above ``threshold`` of the highest bin, with no gap
    longer than ``gap`` between them: any weak patch between the sub-pulses stays below it.

    Parameters
    ----------
    run:
        The run, best that of the masked source, whose profile is smooth.
    bins:
        The bin edges of its map, ``"toa"`` and ``"wavelength"``.
    threshold:
        The fraction of the highest bin a sub-pulse rises above.
    gap:
        A longer gap, in ms, starts the next sub-pulse.

    Returns
    -------
    :
        The time of the highest bin of each sub-pulse, in ms, in order.
    """
    times, [(profile, _)] = _profiles(bins, (run,))
    above = np.flatnonzero(profile > threshold * profile.max())
    groups = np.split(above, np.flatnonzero(np.diff(times[above]) > gap) + 1)
    return [float(times[group[np.argmax(profile[group])]]) for group in groups]


def _set_scale(ax: Axes, log_y: bool, highest: float) -> None:
    """Put the intensity of an axis on a linear scale from 0, or a log one ``LOG_DECADES`` deep.

    Parameters
    ----------
    ax:
        The axis.
    log_y:
        True for the log scale, where the empty bins are left out.
    highest:
        The highest point drawn on it.
    """
    if log_y:  # the empty bins, and the error bars below 0, are left out
        ax.set_yscale("log", nonpositive="mask")
        ax.set_ylim(highest / 10**LOG_DECADES, highest * 2)
    else:
        ax.set_ylim(bottom=0)


def plot_sources(
    bins: Mapping[str, sc.Variable],
    runs: tuple[SourceRun, SourceRun],
    title: str,
    log_y: bool = False,
) -> Figure:
    """Draw the wavelength spectrum and the time profile of both sources.

    Each is summed over the other axis, with its error bars and the points joined by lines.

    Parameters
    ----------
    bins:
        The bin edges of both maps, ``"toa"`` and ``"wavelength"``, evenly spaced.
    runs:
        The runs of both sources, plain then masked.
    title:
        The title of the figure.
    log_y:
        True puts the intensity on a log axis, ``LOG_DECADES`` deep, where the empty bins are
        left out. False, a linear one from 0.

    Returns
    -------
    :
        The figure.
    """
    times, profiles = _profiles(bins, runs)
    wavelengths = sc.midpoints(bins["wavelength"]).to(unit="Å").values
    time_bin, wavelength_bin = times[1] - times[0], wavelengths[1] - wavelengths[0]

    fig, (spectrum_ax, profile_ax) = plt.subplots(1, 2, figsize=(14, 5.5), layout="constrained")
    fig.suptitle(title)
    for run, (profile, profile_error) in zip(runs, profiles, strict=True):
        label = f"{run.component}, error {run.relative_error:.2%}"
        spectrum = run.intensity.sum(axis=0)
        spectrum_error = np.sqrt((run.error**2).sum(axis=0))
        for ax, centers, values, errors in (
            (spectrum_ax, wavelengths, spectrum, spectrum_error),
            (profile_ax, times, profile, profile_error),
        ):
            ax.errorbar(
                centers,
                values,
                yerr=errors,
                label=label,
                fmt=".-",
                markersize=3,
                linewidth=0.8,
                elinewidth=0.8,
            )
    spectrum_ax.set(xlabel="Wavelength [Å]", ylabel=f"Intensity [n/s per {wavelength_bin:.2g} Å]")
    profile_ax.set(xlabel="Time [ms]", ylabel=f"Intensity [n/s per {time_bin:.2g} ms]")
    for ax, highest in (
        (spectrum_ax, max(run.intensity.sum(axis=0).max() for run in runs)),
        (profile_ax, max(profile.max() for profile, _ in profiles)),
    ):
        _set_scale(ax, log_y, highest)
        ax.legend()
    return fig


def plot_peaks(
    bins: Mapping[str, sc.Variable],
    runs: tuple[SourceRun, SourceRun],
    title: str,
    peaks: list[float],
    half_width: float = PEAK_HALF_WIDTH,
) -> Figure:
    """Draw the time profile of both sources around each sub-pulse, a panel each.

    Parameters
    ----------
    bins:
        The bin edges of both maps, ``"toa"`` and ``"wavelength"``, evenly spaced.
    runs:
        The runs of both sources, plain then masked.
    title:
        The title of the figure.
    peaks:
        The times of the sub-pulses, from :func:`find_peaks`, in ms.
    half_width:
        Each panel shows this much on each side of its peak, in ms.

    Returns
    -------
    :
        The figure, ``PEAK_COLUMNS`` panels to a row, in the order of time.
    """
    times, profiles = _profiles(bins, runs)
    time_bin = times[1] - times[0]
    rows = math.ceil(len(peaks) / PEAK_COLUMNS)
    fig, axes = plt.subplots(
        rows,
        PEAK_COLUMNS,
        figsize=(3 * PEAK_COLUMNS, 2.6 * rows),
        layout="constrained",
        squeeze=False,
    )
    fig.suptitle(title)
    for ax, peak in zip(axes.flat, peaks, strict=False):
        inside = np.abs(times - peak) <= half_width
        for run, (profile, error) in zip(runs, profiles, strict=True):
            ax.errorbar(
                times[inside],
                profile[inside],
                yerr=error[inside],
                label=run.component,
                fmt=".-",
                markersize=3,
                linewidth=0.8,
                elinewidth=0.8,
            )
        ax.set(title=f"{peak:.2f} ms", xlim=(peak - half_width, peak + half_width))
        ax.tick_params(labelsize=7)
        ax.set_ylim(bottom=0)
    for ax in axes.flat[len(peaks) :]:  # the empty panels of the last row
        ax.set_axis_off()
    for ax in axes[:, 0]:
        ax.set_ylabel(f"n/s per {time_bin:.2g} ms")
    for ax in axes[-1, :]:
        ax.set_xlabel("Time [ms]")
    axes.flat[0].legend(fontsize=7)
    return fig


def main(
    out: Path | None = None,
    ncount: float = 1e8,
    seed: int = 1,
    lambda_resolution: float = LAMBDA_RESOLUTION,
    time_resolution: float = TIME_RESOLUTION,
) -> list[Path]:
    """Run both sources of :class:`.TRexToSample`, print the table and save the files.

    The figures are ``trex_masked_source.png`` and ``trex_masked_source_log.png``, on a linear
    and on a log intensity axis, from :func:`plot_sources`, and ``trex_masked_source_peaks.png``,
    from :func:`plot_peaks`; ``trex_masked_source.txt`` is the table of :func:`summary_table`.

    Parameters
    ----------
    out:
        Folder of the files. If ``None``, the folder of the plain run, next to its
        ``mccode.h5``.
    ncount:
        Number of McStas rays of each run.
    seed:
        Seed of the random numbers of both runs.
    lambda_resolution:
        The wavelength bin of the sample monitor, in Å.
    time_resolution:
        Its time bin, in s.

    Returns
    -------
    :
        The paths of the figures, linear, log then the peaks, and of the table.
    """
    name = "trex_masked_source"
    plain_instrument = TRexToSample()  # the default control parameters
    masked_instrument = TRexToSample(source_masked=True)
    masked_instrument.name = "trex_to_sample_masked"  # compiled in a folder of its own
    band = wavelength_band(plain_instrument.params.lambda0)
    for instrument in (plain_instrument, masked_instrument):  # in its place, with finer bins
        instrument.sample_monitor = sample_monitor(
            lambda_limits=band,
            lambda_resolution=lambda_resolution,
            time_resolution=time_resolution,
        )
    plain = run_source(plain_instrument, f"{name}_plain", ncount, seed)
    masked = run_source(masked_instrument, f"{name}_masked", ncount, seed)

    table = summary_table(plain, masked, ncount)
    print(table, end="")
    folder = output_folder(out, plain.folder)

    bins = monitor_bins(plain_instrument.sample_monitor)
    peaks = find_peaks(masked, bins)
    conditions = f"ncount = {int(ncount):,} each\n{plain_instrument.running_conditions()}"
    files = [
        save_figure(
            plot_sources(
                bins,
                (plain, masked),
                f"T-Rex flux at the sample monitor, {conditions}",
                log_y=log_y,
            ),
            folder / f"{name}{suffix}.png",
        )
        for suffix, log_y in (("", False), ("_log", True))
    ]
    fig = plot_peaks(
        bins,
        (plain, masked),
        f"The {len(peaks)} sub-pulses at the sample monitor, {conditions}",
        peaks,
    )
    files.append(save_figure(fig, folder / f"{name}_peaks.png"))
    files.append(save_text(table, folder / f"{name}.txt"))
    return files


if __name__ == "__main__":
    main(out=parse_out(__doc__))
