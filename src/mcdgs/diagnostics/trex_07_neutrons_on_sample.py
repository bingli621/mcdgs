"""The neutrons that reach the sample of T-Rex, from McStas and from tof.

:class:`TRexToSample` is :class:`.TRex` up to its sample monitor, at ``sample_pos``: without a
sample, and without the beamstop and the banana detector after it. Its control parameters are
those of :class:`.TRex`.

tof has no guides: every neutron its choppers let through reaches the sample. As in
:mod:`.trex_04_beam_monitor`, its counts are scaled at the source monitor, to the total of
McStas there, where both see the same beam. At the sample monitor, a transmission on top of that
scale stands for what the guides take out of the beam. It is tuned so that tof matches McStas:
by default, so that both have the same total. The transmission is then that of the guides in
McStas, and the shape in time and wavelength that of each model.

The choppers cut the beam into sub-pulses, a bin or two wide and about 4 ms apart.
:func:`.plot_comparison` draws both models twice, with the points of the wavelength spectrum and
the time profile joined by lines: first all on linear scales, then all on log scales, the maps
``LOG_DECADES`` deep, so that the weak bins show too. :func:`find_patches` finds the patches of
each model, the sub-pulses and the weak ones between them, and :func:`plot_patches` draws the
center of mass of each in time and wavelength, with its width as error bars.

:func:`main` runs both models, reads the intensity of each on the bins of the monitors, scales
tof at the source monitor, tunes the transmission at the sample monitor and saves the figures.
Run the module to do it::

    pixi run python -m mcdgs.diagnostics.trex_07_neutrons_on_sample
"""

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import mcstasscript as ms
import numpy as np
import scipp as sc
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from numpy.typing import NDArray

from mcdgs.diagnostics.cli import output_folder, parse_out
from mcdgs.diagnostics.compare import MonitorComparison, plot_comparison, read_comparison
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.mcstas.simulation import Simulation
from mcdgs.tof.builder import build_tof
from mcdgs.tof.reader import monitor_bins

LOG_DECADES = 3  # the log color scale of the maps spans this many decades below the highest bin
PATCH_GAP = 1.0  # ms: a longer gap in time starts a new patch; they are at least 2 ms apart


# --8<-- [start:trex_to_sample]
class TRexToSample(TRex):
    """T-Rex from the source to the sample monitor, without a sample.

    The beamstop and the banana detector, after the sample, are left out.

    Parameters
    ----------
    params:
        Control parameters.
    """

    name = "trex_to_sample"  # the McStas instrument name

    def __init__(self, params: Params = DEFAULT_PARAMS) -> None:
        super().__init__(params=params)
        # after the sample monitor
        del self.beamstop
        del self.banana_detector
        # --8<-- [end:trex_to_sample]


@dataclass(frozen=True)
class Patches:
    """The center of mass of each patch of a time and wavelength map, and its width.

    Each field has one value per patch, in time order. A width is the RMS of the patch about its
    center, with the spread within a bin, its width over √12, added in quadrature: a patch one
    bin wide is that wide, not a point.

    Parameters
    ----------
    time, time_width:
        The mean time of arrival and its width, in ms.
    wavelength, wavelength_width:
        The mean wavelength and its width, in Å.
    intensity:
        The intensity of the patch, in n/s.
    """

    time: NDArray
    time_width: NDArray
    wavelength: NDArray
    wavelength_width: NDArray
    intensity: NDArray


def _mean_and_width(weights: NDArray, edges: NDArray) -> tuple[float, float]:
    """Return the weighted mean of the bin centers and the width, with the spread within a bin."""
    centers = (edges[:-1] + edges[1:]) / 2
    mean = float(np.average(centers, weights=weights))
    variance = np.average((centers - mean) ** 2, weights=weights) + np.diff(edges)[0] ** 2 / 12
    return mean, float(np.sqrt(variance))


def find_patches(
    intensity: NDArray, bins: Mapping[str, sc.Variable], gap: float = PATCH_GAP
) -> Patches:
    """Find the patches of a time and wavelength map: the runs of time bins with neutrons.

    Parameters
    ----------
    intensity:
        The map, with the shape (time bins, wavelength bins).
    bins:
        Its bin edges, ``"toa"`` and ``"wavelength"``, evenly spaced.
    gap:
        The longest gap in time within a patch, in ms; a longer one starts a new patch.

    Returns
    -------
    :
        The center of mass of each patch, its width and its intensity.

    Examples
    --------
    >>> bins = {
    ...     "toa": sc.linspace("toa", 0.0, 6.0, 7, unit="ms"),
    ...     "wavelength": sc.linspace("wavelength", 1.0, 3.0, 3, unit="Å"),
    ... }
    >>> intensity = np.array([[1, 0], [1, 0], [0, 0], [0, 0], [0, 2], [0, 0]])
    >>> patches = find_patches(intensity, bins, gap=1.0)
    >>> patches.time.tolist(), patches.wavelength.tolist(), patches.intensity.tolist()
    ([1.0, 4.5], [1.5, 2.5], [2.0, 2.0])
    """
    times = bins["toa"].to(unit="ms").values
    wavelengths = bins["wavelength"].to(unit="Å").values
    filled = np.flatnonzero(intensity.sum(axis=1) > 0)
    new_patch = np.flatnonzero(np.diff(times[filled]) > gap + np.diff(times)[0]) + 1
    found = []
    for run in np.split(filled, new_patch) if filled.size else []:
        patch = intensity[run[0] : run[-1] + 1]
        time = _mean_and_width(patch.sum(axis=1), times[run[0] : run[-1] + 2])
        wavelength = _mean_and_width(patch.sum(axis=0), wavelengths)
        found.append((*time, *wavelength, float(patch.sum())))
    columns = np.array(found, dtype=float).reshape(-1, 5).T
    return Patches(*columns)


def plot_patches(comparison: MonitorComparison, title: str) -> Figure:
    """Draw the center of mass of each patch of both models, in time and wavelength.

    Each patch is a point at its mean time and wavelength, with its width in each as gray error
    bars, from :func:`find_patches`. Its color is its intensity, on one log scale for both models:
    McStas has small filled circles, tof larger hollow squares with the color as their outline, so
    that where both have a patch the square frames the circle.

    Parameters
    ----------
    comparison:
        Both maps.
    title:
        The title of the figure.

    Returns
    -------
    :
        The figure.
    """
    mcstas = find_patches(comparison.mcstas_intensity, comparison.bins)
    tof = find_patches(comparison.tof_intensity, comparison.bins)
    norm = LogNorm(
        vmin=min(mcstas.intensity.min(), tof.intensity.min()),
        vmax=max(mcstas.intensity.max(), tof.intensity.max()),
    )
    cmap = plt.get_cmap("turbo")

    fig, ax = plt.subplots(figsize=(11, 6), layout="constrained")
    for p in (mcstas, tof):
        ax.errorbar(
            p.time,
            p.wavelength,
            xerr=p.time_width,
            yerr=p.wavelength_width,
            fmt="none",
            ecolor="0.6",
            capsize=2,
            zorder=1,
        )
    # McStas: small filled circles; tof: larger hollow squares, the color as the outline
    points = ax.scatter(
        mcstas.time,
        mcstas.wavelength,
        c=mcstas.intensity,
        cmap=cmap,
        norm=norm,
        marker="o",
        s=25,
        edgecolors="black",
        linewidths=0.3,
        zorder=2,
    )
    ax.scatter(
        tof.time,
        tof.wavelength,
        marker="s",
        s=110,
        facecolors="none",
        edgecolors=cmap(norm(tof.intensity)),
        linewidths=1.8,
        zorder=3,
    )
    # the legend shows the markers only, not the color of one patch
    handles = [
        Line2D(
            [],
            [],
            linestyle="none",
            marker="o",
            markersize=5,
            markerfacecolor="white",
            markeredgecolor="black",
            label=f"McStas, {mcstas.time.size} patches",
        ),
        Line2D(
            [],
            [],
            linestyle="none",
            marker="s",
            markersize=10,
            markerfacecolor="none",
            markeredgecolor="black",
            markeredgewidth=1.8,
            label=f"tof, {tof.time.size} patches",
        ),
    ]
    fig.colorbar(points, ax=ax, label="Intensity of the patch [n/s]")
    ax.set(title=title, xlabel="Time [ms]", ylabel="Wavelength [Å]")
    ax.grid(True, alpha=0.4)
    ax.legend(handles=handles)
    return fig


def main(
    out: Path | None = None,
    ncount: float = 1e8,
    neutrons: int = 10_000_000,
    seed: int = 1,
    transmission: float | None = None,
) -> list[Path]:
    """Run both models of :class:`TRexToSample`, print the transmission and save the figures.

    The figures are ``trex_neutrons_on_sample.png``, on linear scales, and
    ``trex_neutrons_on_sample_log.png``, on log scales, from :func:`.plot_comparison`; and
    ``trex_neutrons_on_sample_patches.png``, the center of mass of each patch, from
    :func:`plot_patches`.

    Parameters
    ----------
    out:
        Folder of the figures. If ``None``, the folder of the McStas run, next to its
        ``mccode.h5``.
    ncount:
        Number of McStas neutrons. Many: the guides and the choppers let through a small
        fraction of the beam.
    neutrons:
        Number of tof neutrons.
    seed:
        Seed of the random numbers of both models.
    transmission:
        The factor on tof at the sample monitor, on top of the scale of the source monitor: the
        fraction of the beam the guides let through, which tof does not model. If ``None``,
        tuned so that tof has the total of McStas at the sample monitor.

    Returns
    -------
    :
        The paths of the figures: linear, log, then the patches.
    """
    instrument = TRexToSample()  # the default control parameters
    name = "trex_neutrons_on_sample"

    # run McStas and load its monitors, in n/s per bin
    folder = Simulation(instrument).run(name, ncount=ncount, seed=seed).parent
    mcstas_data = ms.load_data(str(folder))

    # run tof, then read the monitors of both on the bins of the McStas maps: time, wavelength
    tof_result = build_tof(instrument, neutrons=neutrons, seed=seed).run()
    source_bins = monitor_bins(instrument.source_monitor)
    bins = monitor_bins(instrument.sample_monitor)
    # tof's counts are scaled to the total of McStas at the source monitor, in n/s
    at_source = read_comparison(
        "source_monitor", mcstas_data, tof_result, source_bins, ncount, neutrons
    )

    # at the sample monitor, the same scale, and the transmission on top of it
    comparison = read_comparison(
        "sample_monitor", mcstas_data, tof_result, bins, ncount, neutrons, scale=at_source.scale
    )
    if transmission is None:  # tuned: the transmission that gives tof the total of McStas
        transmission = float(comparison.mcstas_intensity.sum() / comparison.tof_intensity.sum())
    comparison = read_comparison(
        "sample_monitor",
        mcstas_data,
        tof_result,
        bins,
        ncount,
        neutrons,
        scale=at_source.scale,
        transmission=transmission,
    )

    # print the totals, and the scale and transmission of tof
    mcstas_total = comparison.mcstas_intensity.sum()
    tof_counts = comparison.tof_intensity.sum() / comparison.scale / transmission
    print(
        f"McStas total {mcstas_total:.4g} n/s, "
        f"{mcstas_total / at_source.mcstas_intensity.sum():.4g} of that at the source; "
        f"tof {tof_counts:.0f} of {neutrons:,} neutrons, "
        f"{comparison.scale:.4g} n/s each, transmission {transmission:.4g}"
    )

    title = f"T-Rex flux at the sample monitor\n{instrument.running_conditions()}"
    figures = []
    for suffix, log_decades, log_y in (("", None, False), ("_log", LOG_DECADES, True)):
        figure = output_folder(out, folder) / f"{name}{suffix}.png"
        fig = plot_comparison(comparison, title, log_decades=log_decades, connect=True, log_y=log_y)
        fig.savefig(figure, dpi=150)
        plt.close(fig)
        print(f"Saved {figure}")
        figures.append(figure)

    figure = output_folder(out, folder) / f"{name}_patches.png"
    fig = plot_patches(
        comparison, f"Center of mass of each patch\n{instrument.running_conditions()}"
    )
    fig.savefig(figure, dpi=150)
    plt.close(fig)
    print(f"Saved {figure}")
    figures.append(figure)
    return figures


if __name__ == "__main__":
    main(out=parse_out(__doc__))
