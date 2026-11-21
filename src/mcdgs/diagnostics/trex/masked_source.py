"""The neutrons that reach the sample of T-Rex, from ESS_butterfly and from Polygon_ESS_butterfly.

Both runs are :class:`TRexToSampleMonitor`, T-Rex up to its sample monitor, with the same number
of rays and the same seed. The first
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
``PEAK_HALF_WIDTH`` on each side, both with :func:`.plot_profile`. :func:`summary_table` lists, for each source, the total at
the sample monitor and its relative error, the time of the run and the figure of merit, one over
the relative error squared times the time: how much a second of each run is worth.

:func:`main` runs both and saves the figures and the table. Run the module to do it::

    pixi run python -m mcdgs.diagnostics.trex.masked_source
"""

import json
import math
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipp as sc
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from numpy.typing import NDArray

from mcdgs.diagnostics import cli, compare
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.instruments.trex.monitors import sample_monitor
from mcdgs.instruments.trex.sources import wavelength_band
from mcdgs.mcstas import Simulation, read_mcstas_map
from mcdgs.reduction.binning import monitor_bins

MONITOR = "sample_monitor"  # where both sources are compared
LAMBDA_RESOLUTION = 0.002  # Å, the wavelength bin of its map: a tenth of that of T-Rex
TIME_RESOLUTION = 1e-5  # s, its time bin: a fifth of that of T-Rex
PEAK_THRESHOLD = 0.01  # a sub-pulse has bins above this fraction of the highest bin of the profile
PEAK_GAP = 1.0  # ms, a gap longer than this between such bins starts the next sub-pulse
PEAK_HALF_WIDTH = 0.1  # ms, the zoom on a sub-pulse shows this much on each side of its peak
PEAK_COLUMNS = 6  # the zooms on the sub-pulses are drawn this many to a row
LOG_DECADES = 6  # a log axis spans this many decades below the highest point


# --8<-- [start:instrument]
class TRexToSampleMonitor(TRex):
    """T-Rex from the source to its sample monitor, without a sample, with finer bins there.

    The beamstop and the banana detector, after the sample, are left out. The sample monitor
    spans the band of the source, with time and wavelength bins finer than those of T-Rex, to
    show the shape of each sub-pulse.

    Parameters
    ----------
    params:
        Control parameters.
    source_masked:
        As in :class:`.TRex`: False for the plain source, True for the masked one.
    lambda_resolution:
        The wavelength bin of the sample monitor, in Å.
    time_resolution:
        Its time bin, in s.
    """

    name = "trex_to_sample_monitor"  # the McStas instrument name

    def __init__(
        self,
        params: Params = DEFAULT_PARAMS,
        source_masked: bool = False,
        lambda_resolution: float = LAMBDA_RESOLUTION,
        time_resolution: float = TIME_RESOLUTION,
    ) -> None:
        super().__init__(params=params, source_masked=source_masked)
        # after the sample monitor
        del self.beamstop
        del self.banana_detector
        self.sample_monitor = sample_monitor(  # in its place, with finer bins
            lambda_limits=wavelength_band(params.lambda0),
            lambda_resolution=lambda_resolution,
            time_resolution=time_resolution,
        )


def define_instrument(
    source_masked: bool = False,
    lambda_resolution: float = LAMBDA_RESOLUTION,
    time_resolution: float = TIME_RESOLUTION,
) -> TRexToSampleMonitor:
    """Return an instrument of this diagnostic: :class:`TRexToSampleMonitor` with one source.

    The masked one has a name of its own, so that each is compiled in a build folder of its own.

    Parameters
    ----------
    source_masked:
        True for the masked source, False for the plain one.
    lambda_resolution, time_resolution:
        As in :class:`TRexToSampleMonitor`.

    Returns
    -------
    :
        The instrument.
    """
    params = Params()  # the default control parameters
    instrument = TRexToSampleMonitor(
        params=params,
        source_masked=source_masked,
        lambda_resolution=lambda_resolution,
        time_resolution=time_resolution,
    )
    if source_masked:
        instrument.name = "trex_to_sample_monitor_masked"
    return instrument


# --8<-- [end:instrument]


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
    sample_map:
        The map of the sample monitor, as :func:`.read_mcstas_map`: in n/s per bin, (time bins,
        wavelength bins), with its bin edges, ``"toa"`` and ``"wavelength"``, evenly spaced, and
        its errors as variances.
    acceptance:
        The acceptance of the masked source, from its region file, or ``None`` for the plain one.
    """

    component: str
    folder: Path
    seconds: float
    sample_map: sc.DataArray
    acceptance: float | None

    @property
    def total(self) -> float:
        """The total intensity at the sample monitor, in n/s."""
        return float(self.sample_map.sum().value)

    @property
    def relative_error(self) -> float:
        """The relative error of :attr:`total`: the errors of the bins add in quadrature."""
        return float(sc.stddevs(self.sample_map.sum()).value) / self.total

    @property
    def spectrum(self) -> sc.DataArray:
        """The wavelength spectrum: the map summed over time, its errors in quadrature."""
        return self.sample_map.sum("toa")

    @property
    def time_profile(self) -> sc.DataArray:
        """The time profile: the map summed over wavelength, its errors in quadrature."""
        return self.sample_map.sum("wavelength")

    @property
    def figure_of_merit(self) -> float:
        """One over the relative error squared times the time, in 1/s."""
        return 1 / (self.relative_error**2 * self.seconds)


def run_source(
    instrument: TRexToSampleMonitor, folder_name: str, ncount: float, seed: int
) -> SourceRun:
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
    sim = Simulation(instrument)
    run = sim.run(folder_name=folder_name, ncount=ncount, seed=seed)
    seconds = time.perf_counter() - start
    sample_map = read_mcstas_map(MONITOR, run, monitor_bins(instrument.sample_monitor.bin_limits))
    region = run.folder / f"{instrument.source.name}.json"  # the masked source writes its region
    acceptance = json.loads(region.read_text())["acceptance"] if region.is_file() else None
    return SourceRun(instrument.source.component_name, run.folder, seconds, sample_map, acceptance)


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
    return cli.markdown_table(header, rows)


def _times(profile: sc.DataArray) -> NDArray:
    """Return the centers of the time bins of a time profile, in ms."""
    return sc.midpoints(profile.coords["toa"]).to(unit="ms").values


def find_peaks(
    run: SourceRun,
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
    threshold:
        The fraction of the highest bin a sub-pulse rises above.
    gap:
        A longer gap, in ms, starts the next sub-pulse.

    Returns
    -------
    :
        The time of the highest bin of each sub-pulse, in ms, in order.
    """
    times, profile = _times(run.time_profile), run.time_profile.values
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


def plot_sources(runs: tuple[SourceRun, SourceRun], title: str, log_y: bool = False) -> Figure:
    """Draw the wavelength spectrum and the time profile of both sources.

    Each is summed over the other axis, with its error bars and the points joined by lines, as
    :func:`.plot_profile` draws them.

    Parameters
    ----------
    runs:
        The runs of both sources, plain then masked, on the same bins.
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
    fig, (spectrum_ax, profile_ax) = plt.subplots(1, 2, figsize=(14, 5.5), layout="constrained")
    fig.suptitle(title)
    for run in runs:
        label = f"{run.component}, error {run.relative_error:.2%}"
        compare.plot_profile(spectrum_ax, run.spectrum, label, connect=True)
        compare.plot_profile(profile_ax, run.time_profile, label, connect=True)
    for ax, highest in (
        (spectrum_ax, max(run.spectrum.values.max() for run in runs)),
        (profile_ax, max(run.time_profile.values.max() for run in runs)),
    ):
        _set_scale(ax, log_y, highest)
        ax.legend()
    return fig


def plot_peaks(
    runs: tuple[SourceRun, SourceRun],
    title: str,
    peaks: list[float],
    half_width: float = PEAK_HALF_WIDTH,
) -> Figure:
    """Draw the time profile of both sources around each sub-pulse, a panel each.

    Each panel has the time bins within ``half_width`` of its peak, as :func:`.plot_profile`
    draws them, with their axis labels on the outer panels only.

    Parameters
    ----------
    runs:
        The runs of both sources, plain then masked, on the same bins.
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
    times = _times(runs[0].time_profile)
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
        inside = np.flatnonzero(np.abs(times - peak) <= half_width)  # the bins around the peak
        for run in runs:
            around = run.time_profile["toa", inside[0] : inside[-1] + 1]
            compare.plot_profile(ax, around, run.component, connect=True)
        ax.set(title=f"{peak:.2f} ms", xlim=(peak - half_width, peak + half_width))
        ax.tick_params(labelsize=7)
        ax.set_ylim(bottom=0)
    for ax in axes[:, 1:].flat:  # the intensity is labeled on the first column only
        ax.set_ylabel("")
    for ax in axes[:-1, :].flat:  # and the time on the last row only
        ax.set_xlabel("")
    for ax in axes.flat[len(peaks) :]:  # the empty panels of the last row
        ax.set_axis_off()
    axes.flat[0].legend(fontsize=7)
    return fig


def main(
    out: Path | None = None,
    ncount: float = 1e8,
    seed: int = 1,
    lambda_resolution: float = LAMBDA_RESOLUTION,
    time_resolution: float = TIME_RESOLUTION,
) -> list[Path]:
    """Run both sources of :class:`TRexToSampleMonitor`, print the table and save the files.

    The figures are ``trex_masked_source.png`` and ``trex_masked_source_log.png``, on a linear
    and on a log intensity axis, from :func:`plot_sources`, and ``trex_masked_source_peaks.png``,
    from :func:`plot_peaks`; ``trex_masked_source.txt`` is the table of :func:`summary_table`.

    Parameters
    ----------
    out:
        Folder of the files. If ``None``, ``trex_masked_source`` in the results folder of
        :meth:`.Workspace.find`.
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
    plain_instrument, masked_instrument = (
        define_instrument(masked, lambda_resolution, time_resolution) for masked in (False, True)
    )
    plain = run_source(plain_instrument, f"{name}_plain", ncount, seed)
    masked = run_source(masked_instrument, f"{name}_masked", ncount, seed)

    table = summary_table(plain, masked, ncount)
    print(table, end="")
    folder = cli.output_folder(out, "trex_masked_source")

    peaks = find_peaks(masked)
    conditions = f"ncount = {int(ncount):,} each\n{plain_instrument.running_conditions()}"
    files = [
        cli.save_figure(
            plot_sources(
                (plain, masked),
                f"T-Rex flux at the sample monitor, {conditions}",
                log_y=log_y,
            ),
            folder / f"{name}{suffix}.png",
        )
        for suffix, log_y in (("", False), ("_log", True))
    ]
    fig = plot_peaks(
        (plain, masked),
        f"The {len(peaks)} sub-pulses at the sample monitor, {conditions}",
        peaks,
    )
    files.append(cli.save_figure(fig, folder / f"{name}_peaks.png"))
    files.append(cli.save_text(table, folder / f"{name}.txt"))
    return files


if __name__ == "__main__":
    main(out=cli.parse_out(__doc__))
