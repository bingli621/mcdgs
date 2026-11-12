"""The neutrons that reach the sample of T-Rex, from McStas and from tof.

:class:`TRexToSample` is :class:`.TRex` up to its sample monitor, at ``sample_pos``: without a
sample, and without the beamstop and the banana detector after it. Its control parameters are
those of :class:`.TRex`. After the sample monitor it has :func:`sample_xy_monitor`, in the same
place, which records where the neutrons cross the sample position: one map of x against y.
:func:`read_xy_map` reads that map from the McStas run and :func:`plot_xy_map` draws it, with
its profile along each axis. tof has no x and y, so the map is that of McStas alone.

tof has no guides: every neutron its choppers let through reaches the sample. As in
:mod:`.diagnostics.trex.beam_monitor`, its counts are scaled at the source monitor, to the total of
McStas there, where both see the same beam. At the sample monitor, a transmission on top of that
scale stands for what the guides take out of the beam. It is tuned so that tof matches McStas:
by default, so that both have the same total. The transmission is then that of the guides in
McStas, and the shape in time and wavelength that of each model.

The choppers cut the beam into sub-pulses, a bin or two wide and about 4 ms apart.
:func:`.plot_comparison` draws both models twice, with the points of the wavelength spectrum and
the time profile joined by lines: first all on linear scales, then all on log scales, the maps
``LOG_DECADES`` deep, so that the weak bins show too. :func:`.find_patches` finds the patches of
each model, the sub-pulses and any weak ones between them, and :func:`plot_patches` draws the
center of mass of each in time and wavelength, with its width as error bars.

:func:`subpulse_windows` reduces the tof run as :mod:`.diagnostics.trex.elastic_plane` reduces McStas:
the sub-pulses at the pre-sample and beamstop monitors give their incident energy and the time
they reach the sample, and :func:`.arrival_window` when the neutrons of each reach the banana
detector. :func:`windows_table` lists them.

:func:`main` runs both models, reads the intensity of each on the bins of the monitors, scales
tof at the source monitor, tunes the transmission at the sample monitor and saves the figures,
that of the x and y map among them.
Run the module to do it::

    pixi run python -m mcdgs.diagnostics.trex.neutrons_on_sample
"""

from pathlib import Path

import matplotlib.pyplot as plt
import mcstasscript as ms
import numpy as np
import scipp as sc
import tof
from matplotlib.artist import Artist
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mcstasscript.data.data import McStasData
from scippneutron.conversion import tof as tof_kernels

from mcdgs.core.components import Monitor_nD
from mcdgs.core.utils import bins_from_limits
from mcdgs.diagnostics.cli import markdown_table, output_folder, parse_out, save_figure, save_text
from mcdgs.diagnostics.compare import MonitorComparison, plot_comparison, read_scaled_comparisons
from mcdgs.instruments.trex.detectors import BANANA_RADIUS
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.instruments.trex.monitors import SAMPLE_MONITOR_HEIGHT, SAMPLE_MONITOR_WIDTH
from mcdgs.mcstas.simulation import Simulation
from mcdgs.reduction.inelastic import DEFAULT_ENERGY_TRANSFER_RATIO, arrival_window, time_at_sample
from mcdgs.reduction.patches import find_patches, subpulse_flight
from mcdgs.tof.builder import build_tof
from mcdgs.tof.reader import monitor_bins, read_tof_monitor

LOG_DECADES = 3  # the log color scale of the maps spans this many decades below the highest bin
MONITORS = ("pre_sample_monitor", "beamstop_monitor")  # give the incident energy, in beam order
SAMPLE_XY_MONITOR_RESOLUTION = 0.002  # m, the width of a bin of x and of y, 2 mm


def sample_xy_monitor(resolution: float = SAMPLE_XY_MONITOR_RESOLUTION) -> Monitor_nD:
    """Create the x and y monitor at the sample position.

    It records where the neutrons cross the plane of ``sample_pos``, across the beam, as one 2D
    map of x against y, summed over time and wavelength, and leaves the neutrons unchanged. It
    has the place and the opening of the sample monitor, ``SAMPLE_MONITOR_WIDTH`` by
    ``SAMPLE_MONITOR_HEIGHT``, centered on the beam axis: x is horizontal, y is up.

    Parameters
    ----------
    resolution:
        Width of a bin of x and of y, in m. Each axis has the fewest bins that cover the
        opening, see :func:`.bins_from_limits`.

    Returns
    -------
    :
        A new monitor.
    """
    x_start, x_stop, x_bins = bins_from_limits(
        (-SAMPLE_MONITOR_WIDTH / 2, SAMPLE_MONITOR_WIDTH / 2), resolution
    )
    y_start, y_stop, y_bins = bins_from_limits(
        (-SAMPLE_MONITOR_HEIGHT / 2, SAMPLE_MONITOR_HEIGHT / 2), resolution
    )
    return Monitor_nD(
        xwidth=SAMPLE_MONITOR_WIDTH,
        yheight=SAMPLE_MONITOR_HEIGHT,
        restore_neutron=1,
        options=(  # x and y in m: one 2D map
            f"x limits=[{x_start} {x_stop}] bins={x_bins} "
            f"y limits=[{y_start} {y_stop}] bins={y_bins}"
        ),
        filename="sample_xy_monitor",
        relative="sample_pos",
    )


# --8<-- [start:trex_to_sample]
class TRexToSample(TRex):
    """T-Rex from the source to the sample monitor, without a sample.

    The beamstop and the banana detector, after the sample, are left out. In their place is
    :func:`sample_xy_monitor`, at ``sample_pos`` as the sample monitor.

    Parameters
    ----------
    params:
        Control parameters.
    source_masked:
        As in :class:`.TRex`.
    """

    name = "trex_to_sample"  # the McStas instrument name

    def __init__(self, params: Params = DEFAULT_PARAMS, source_masked: bool = False) -> None:
        super().__init__(params=params, source_masked=source_masked)
        # after the sample monitor
        del self.beamstop
        del self.banana_detector
        self.sample_xy_monitor = sample_xy_monitor()
        # --8<-- [end:trex_to_sample]


def read_xy_map(monitor: Monitor_nD, mcstas_data: list[McStasData]) -> sc.DataArray:
    """Read the x and y map of a monitor of one McStas run.

    Parameters
    ----------
    monitor:
        The monitor, such as :func:`sample_xy_monitor`: its ``filename`` is its name in the run,
        and its ``options`` give one 2D map of x against y, in m.
    mcstas_data:
        The monitors of the McStas run, as ``mcstasscript.load_data`` reads them.

    Returns
    -------
    :
        The intensity per bin, in n/s, with dims ``y`` and ``x`` and their bin edges, in m. Its
        variances are the squares of the errors of McStas.

    Raises
    ------
    ValueError
        If the map of McStas does not have the bins of ``options``.
    """
    name = str(monitor.filename)
    data = ms.name_search(name, mcstas_data)
    limits = monitor.bin_limits
    edges = {
        axis: sc.linspace(axis, start, stop, n_bins + 1, unit="m")
        for axis, (start, stop, n_bins) in ((axis, limits[axis]) for axis in ("y", "x"))
    }
    # McStasScript has the map as (y, x), the y axis first
    intensity = np.asarray(data.Intensity, dtype=float)
    shape = tuple(len(edge) - 1 for edge in edges.values())
    if intensity.shape != shape:
        raise ValueError(
            f"{name}: McStas has a map of {intensity.shape} bins, the bins are {shape}"
        )
    return sc.DataArray(
        sc.array(
            dims=["y", "x"],
            values=intensity,
            variances=np.asarray(data.Error, dtype=float) ** 2,
            unit="1/s",
        ),
        coords=edges,
    )


def plot_xy_map(xy_map: sc.DataArray, title: str) -> Figure:
    """Draw the x and y map of the beam, with its profile along each axis.

    The map is on the left, in n/s per bin, with both axes to the same scale. On the right are
    its profiles: the map summed over y, against x, and summed over x, against y, each with the
    errors of McStas and, in its legend, its full width at half maximum.

    Parameters
    ----------
    xy_map:
        The intensity per bin, as :func:`read_xy_map`.
    title:
        The title of the figure.

    Returns
    -------
    :
        The figure.
    """
    fig = plt.figure(figsize=(12, 6), layout="constrained")
    fig.suptitle(title)
    map_ax, x_ax, y_ax = fig.subplot_mosaic([["map", "x"], ["map", "y"]]).values()

    mesh = map_ax.pcolormesh(
        xy_map.coords["x"].to(unit="mm").values,
        xy_map.coords["y"].to(unit="mm").values,
        xy_map.values,
        cmap="turbo",
        vmin=0,
    )
    map_ax.set(xlabel="x [mm], horizontal", ylabel="y [mm], up", aspect="equal")
    fig.colorbar(mesh, ax=map_ax, label="Intensity [n/s per bin]")

    for ax, axis, summed, label in ((x_ax, "x", "y", "x [mm]"), (y_ax, "y", "x", "y [mm]")):
        profile = xy_map.sum(summed)
        centers = sc.midpoints(profile.coords[axis]).to(unit="mm").values
        ax.errorbar(
            centers,
            profile.values,
            yerr=np.sqrt(profile.variances),
            drawstyle="steps-mid",
            linewidth=1,
            label=f"FWHM {full_width_at_half_maximum(centers, profile.values):.1f} mm",
        )
        ax.set(title=f"Summed over {summed}", xlabel=label, ylabel="Intensity [n/s per bin]")
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.4)
        ax.legend()
    return fig


def full_width_at_half_maximum(centers: np.ndarray, values: np.ndarray) -> float:
    """Return the full width at half maximum of a profile with one peak.

    The width runs between the outermost crossings of half the highest value, each found on the
    straight line between the two points around it.

    Parameters
    ----------
    centers:
        The centers of the bins, rising.
    values:
        The value of each bin. The first and the last must be below half the highest.

    Returns
    -------
    :
        The width, in the unit of ``centers``.

    Raises
    ------
    ValueError
        If the first or the last value is not below half the highest: the peak is cut off.

    Examples
    --------
    >>> full_width_at_half_maximum(np.array([0.0, 1.0, 2.0, 3.0, 4.0]), np.array([0, 2, 4, 2, 0]))
    2.0
    """
    half = values.max() / 2
    above = np.flatnonzero(values >= half)
    first, last = above[0], above[-1]
    if first == 0 or last == len(values) - 1:
        raise ValueError("the profile does not fall below half its maximum at both ends")
    left = np.interp(half, values[[first - 1, first]], centers[[first - 1, first]])
    right = np.interp(half, values[[last + 1, last]], centers[[last + 1, last]])
    return float(right - left)


def plot_patches(
    comparison: MonitorComparison, title: str, windows: sc.DataArray | None = None
) -> Figure:
    """Draw the center of mass of each patch of both models, in time and wavelength.

    Each patch is a point at its mean time and wavelength, with its width in each as gray error
    bars, from :func:`.find_patches`. Its color is its intensity, on one log scale for both models:
    McStas has small filled circles, tof larger hollow squares with the color as their outline, so
    that where both have a patch the square frames the circle.

    With ``windows``, each sub-pulse is also drawn at its incident wavelength: a black cross when
    it reaches the sample, and its arrival window at the banana detector, from the opening to the
    close, as a gray band across all wavelengths, with a thin line when its neutrons arrive if
    they keep their speed through the sample.

    Parameters
    ----------
    comparison:
        Both maps.
    title:
        The title of the figure.
    windows:
        The sub-pulses of :func:`subpulse_windows`, or ``None``.

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
    handles: list[Artist] = [
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
    if windows is not None:
        wavelength = windows.coords["incident_wavelength"].to(unit="Å").values
        ax.scatter(
            windows.coords["time_at_sample"].to(unit="ms").values,
            wavelength,
            marker="x",
            s=40,
            color="black",
            linewidths=1.2,
            zorder=4,
        )
        for opens, closes in zip(
            windows.coords["earliest_arrival"].to(unit="ms").values,
            windows.coords["latest_arrival"].to(unit="ms").values,
        ):
            ax.axvspan(opens, closes, color="0.5", alpha=0.2, linewidth=0, zorder=0)
        for elastic in windows.coords["elastic_arrival"].to(unit="ms").values:
            ax.axvline(elastic, color="0.3", linewidth=0.6, zorder=0.5)
        handles += [
            Line2D([], [], linestyle="none", marker="x", color="black", label="At the sample"),
            Patch(
                facecolor="0.5",
                alpha=0.2,
                label="Window at the detector, E/Ei from \u22120.8 to 0.8",
            ),
            Line2D([], [], color="0.3", linewidth=0.6, label="Elastic arrival at the detector"),
        ]
    fig.colorbar(points, ax=ax, label="Intensity of the patch [n/s]")
    ax.set(title=title, xlabel="Time [ms]", ylabel="Wavelength [Å]")
    ax.grid(True, alpha=0.4)
    ax.legend(handles=handles)
    return fig


def subpulse_windows(
    instrument: TRex,
    tof_result: tof.Result,
    energy_transfer_ratio: tuple[float, float] = DEFAULT_ENERGY_TRANSFER_RATIO,
) -> sc.DataArray:
    """Give each sub-pulse of a tof run its arrival window at the banana detector.

    As :mod:`.diagnostics.trex.elastic_plane` does with McStas: :func:`.find_patches` finds the
    sub-pulses at both monitors, in the first pulse, :func:`.subpulse_flight` gives their
    incident energy from the flight between them, :func:`.time_at_sample` the time they reach
    the sample, and :func:`.arrival_window` when their neutrons reach the detector,
    ``BANANA_RADIUS`` from the sample.

    Parameters
    ----------
    instrument:
        The instrument of the run, with the monitors of ``MONITORS``.
    tof_result:
        The tof run.
    energy_transfer_ratio:
        The smallest and the largest ``E / Ei`` of the windows, as :func:`.arrival_window`.

    Returns
    -------
    :
        The sub-pulses along ``subpulse``, with the coordinates of :func:`.subpulse_flight`,
        ``time_at_sample``, and ``earliest_arrival``, ``elastic_arrival``, of a neutron that
        keeps its speed, and ``latest_arrival`` at the detector, in ms. The data are tof
        neutrons at the pre-sample monitor, not n/s.
    """
    positions = {name: sc.vector(xyz, unit="m") for name, xyz in instrument.positions().items()}
    start, end = (
        find_patches(read_tof_monitor(tof_result[name], bins), bins)
        for name in MONITORS
        for bins in [monitor_bins(getattr(instrument, name))]
    )
    pulses = subpulse_flight(start, positions[MONITORS[0]], end, positions[MONITORS[1]])
    pulses.unit = None  # tof neutrons, not n/s
    pulses.coords["time_at_sample"] = time_at_sample(
        pulses.coords["time_at_start"],
        positions[MONITORS[0]],
        positions["sample_pos"],
        pulses.coords["incident_energy"],
    )
    window = arrival_window(
        pulses.coords["time_at_sample"],
        sc.scalar(BANANA_RADIUS, unit="m"),
        pulses.coords["incident_energy"],
        energy_transfer_ratio=energy_transfer_ratio,
    )
    pulses.coords.update(window)
    # a neutron that flies through the sample without changing its speed
    at_sample = pulses.coords["time_at_sample"]
    flight = tof_kernels.tof_from_wavelength(
        wavelength=pulses.coords["incident_wavelength"], Ltotal=sc.scalar(BANANA_RADIUS, unit="m")
    )
    pulses.coords["elastic_arrival"] = at_sample + flight.to(unit=at_sample.unit)
    return pulses


def windows_table(windows: sc.DataArray) -> str:
    """Return the arrival windows of :func:`subpulse_windows` as a Markdown table.

    Each row is a sub-pulse: its incident wavelength and energy, when it reaches the sample,
    and when its window at the detector opens, its elastic neutrons arrive and the window
    closes. The last column is the time from the close of its window to the opening of the
    next: negative where the two overlap.

    Parameters
    ----------
    windows:
        The sub-pulses of :func:`subpulse_windows`.

    Returns
    -------
    :
        The table.
    """
    wavelength = windows.coords["incident_wavelength"].to(unit="Å").values
    energy = windows.coords["incident_energy"].to(unit="meV").values
    at_sample = windows.coords["time_at_sample"].to(unit="ms").values
    opens = windows.coords["earliest_arrival"].to(unit="ms").values
    elastic = windows.coords["elastic_arrival"].to(unit="ms").values
    closes = windows.coords["latest_arrival"].to(unit="ms").values
    header = [
        "Sub-pulse",
        "λi",
        "Ei",
        "At the sample",
        "Window opens",
        "Elastic",
        "Window closes",
        "To the next",
    ]
    to_next = [f"{after - close:.3f} ms" for close, after in zip(closes, opens[1:])] + [""]
    rows = [
        (
            i,
            f"{wavelength[i]:.4f} Å",
            f"{energy[i]:.3f} meV",
            f"{at_sample[i]:.3f} ms",
            f"{opens[i]:.3f} ms",
            f"{elastic[i]:.3f} ms",
            f"{closes[i]:.3f} ms",
            to_next[i],
        )
        for i in range(windows.sizes["subpulse"])
    ]
    return markdown_table(header, rows)


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
    :func:`plot_patches`. ``trex_neutrons_on_sample_windows.txt`` is the table of the arrival
    windows of tof, from :func:`windows_table`. ``trex_neutrons_on_sample_xy.png`` is the x and y
    map of McStas at the sample position, from :func:`plot_xy_map`.

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
        The paths of the figures: linear, log, then the patches; of the table of the
        arrival windows; and of the x and y map.
    """
    instrument = TRexToSample()  # the default control parameters
    name = "trex_neutrons_on_sample"

    # run McStas and load its monitors, in n/s per bin
    folder = Simulation(instrument).run(name, ncount=ncount, seed=seed).parent
    mcstas_data = ms.load_data(str(folder))

    # run tof, then read the monitors of both on the bins of the McStas maps: time, wavelength
    tof_result = build_tof(instrument, neutrons=neutrons, seed=seed).run()
    # tof's counts are scaled to the total of McStas at the source monitor, in n/s; at the
    # sample monitor, the same scale, and the transmission on top of it
    comparisons = read_scaled_comparisons(
        mcstas_data,
        tof_result,
        {
            "source_monitor": monitor_bins(instrument.source_monitor),
            "sample_monitor": monitor_bins(instrument.sample_monitor),
        },
        ncount,
        neutrons,
        transmission=transmission,
    )
    at_source, comparison = comparisons["source_monitor"], comparisons["sample_monitor"]
    transmission = comparison.transmission

    # print the totals, and the scale and transmission of tof
    mcstas_total = comparison.mcstas_intensity.sum()
    tof_counts = comparison.tof_intensity.sum() / comparison.scale / transmission
    print(
        f"McStas total {mcstas_total:.4g} n/s, "
        f"{mcstas_total / at_source.mcstas_intensity.sum():.4g} of that at the source; "
        f"tof {tof_counts:.0f} of {neutrons:,} neutrons, "
        f"{comparison.scale:.4g} n/s each, transmission {transmission:.4g}"
    )

    folder = output_folder(out, folder)
    conditions = instrument.running_conditions()
    files = [
        save_figure(
            plot_comparison(
                comparison,
                f"T-Rex flux at the sample monitor\n{conditions}",
                log_decades=log_decades,
                connect=True,
                log_y=log_y,
            ),
            folder / f"{name}{suffix}.png",
        )
        for suffix, log_decades, log_y in (("", None, False), ("_log", LOG_DECADES, True))
    ]

    windows = subpulse_windows(instrument, tof_result)
    fig = plot_patches(comparison, f"Center of mass of each patch\n{conditions}", windows=windows)
    files.append(save_figure(fig, folder / f"{name}_patches.png"))

    table = windows_table(windows)
    print(table, end="")
    files.append(save_text(table, folder / f"{name}_windows.txt"))

    xy_map = read_xy_map(instrument.sample_xy_monitor, mcstas_data)
    fig = plot_xy_map(xy_map, f"T-Rex beam at the sample position, McStas\n{conditions}")
    files.append(save_figure(fig, folder / f"{name}_xy.png"))
    return files


if __name__ == "__main__":
    main(out=parse_out(__doc__))
