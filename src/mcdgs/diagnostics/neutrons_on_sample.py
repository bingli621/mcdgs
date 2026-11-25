"""The neutrons that reach the sample, read and drawn, for any instrument.

:func:`sample_xy_monitor` records where the neutrons cross the sample position, one map of x
against y; :func:`read_xy_map` reads that map from a McStas run and :func:`plot_xy_map` draws it,
with its profile along each axis. :func:`plot_patches` draws the center of mass of each patch at
the sample monitor, in McStas and tof, with the arrival window of each sub-pulse at the detector,
:func:`elastic_arrival` gives when the neutrons of a sub-pulse that keep their speed reach the
detector, and :func:`windows_table` lists the windows. The neutrons on sample diagnostic of each
instrument finds its sub-pulses and their windows.
"""

import matplotlib.pyplot as plt
import numpy as np
import scipp as sc
from matplotlib.artist import Artist
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scippneutron.conversion import tof as tof_kernels

from mcdgs.core import bins_from_limits
from mcdgs.core.components import Monitor_nD
from mcdgs.diagnostics import cli, compare
from mcdgs.mcstas import Run, read_mcstas_maps
from mcdgs.reduction.patches import find_patches

SAMPLE_XY_MONITOR_RESOLUTION = 0.002  # m, the width of a bin of x and of y, 2 mm


def sample_xy_monitor(
    width: float, height: float, resolution: float = SAMPLE_XY_MONITOR_RESOLUTION
) -> Monitor_nD:
    """Create the x and y monitor at the sample position.

    It records where the neutrons cross the plane of ``sample_pos``, across the beam, as one 2D
    map of x against y, summed over time and wavelength, and leaves the neutrons unchanged. It
    is relative to ``sample_pos``, the place of the sample monitor, and centered on the beam
    axis: x is horizontal, y is up.

    Parameters
    ----------
    width, height:
        Size of the opening, in m, that of the sample monitor.
    resolution:
        Width of a bin of x and of y, in m. Each axis has the fewest bins that cover the
        opening, see :func:`.bins_from_limits`.

    Returns
    -------
    :
        A new monitor.
    """
    x_start, x_stop, x_bins = bins_from_limits((-width / 2, width / 2), resolution)
    y_start, y_stop, y_bins = bins_from_limits((-height / 2, height / 2), resolution)
    return Monitor_nD(
        xwidth=width,
        yheight=height,
        restore_neutron=1,
        options=(  # x and y in m: one 2D map
            f"x limits=[{x_start} {x_stop}] bins={x_bins} "
            f"y limits=[{y_start} {y_stop}] bins={y_bins}"
        ),
        filename="sample_xy_monitor",
        relative="sample_pos",
    )


def read_xy_map(monitor: Monitor_nD, run: Run) -> sc.DataArray:
    """Read the x and y map of a monitor of one McStas run.

    Parameters
    ----------
    monitor:
        The monitor, such as :func:`sample_xy_monitor`: its ``filename`` is its name in the run,
        and its ``options`` give one 2D map of x against y, in m.
    run:
        The McStas run.

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
    limits = monitor.bin_limits
    edges = {  # in the order of the axes of options, x then y
        axis: sc.linspace(axis, start, stop, n_bins + 1, unit="m")
        for axis, (start, stop, n_bins) in ((axis, limits[axis]) for axis in ("x", "y"))
    }
    name = str(monitor.filename)
    return read_mcstas_maps(run, {name: edges})[name].transpose(["y", "x"]).copy()


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
    comparison: compare.MonitorComparison, title: str, windows: sc.DataArray | None = None
) -> Figure:
    """Draw the center of mass of each patch of both models, in time and wavelength.

    Each patch is a point at its mean time and wavelength, with its width in each as gray error
    bars, from :func:`.find_patches`. Its color is its intensity, on one log scale for both models:
    McStas has small filled circles, tof larger hollow squares with the color as their outline, so
    that where both have a patch the square frames the circle.

    With ``windows``, each sub-pulse is also drawn at its incident wavelength: a black cross when
    it reaches the sample, and its arrival window at the detector, from the opening to the
    close, as a gray band across all wavelengths, with a thin line when its neutrons arrive if
    they keep their speed through the sample.

    Parameters
    ----------
    comparison:
        Both maps.
    title:
        The title of the figure.
    windows:
        The sub-pulses with their arrival windows, or ``None``: along ``subpulse``, with the
        coordinates of :func:`windows_table`.

    Returns
    -------
    :
        The figure.
    """
    mcstas = find_patches(comparison.first.values, comparison.bins)
    tof = find_patches(comparison.second.values, comparison.bins)
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


def windows_table(windows: sc.DataArray) -> str:
    """Return the arrival windows of the sub-pulses as a Markdown table.

    Each row is a sub-pulse: its incident wavelength and energy, when it reaches the sample,
    and when its window at the detector opens, its elastic neutrons arrive and the window
    closes. The last column is the time from the close of its window to the opening of the
    next: negative where the two overlap.

    Parameters
    ----------
    windows:
        The sub-pulses along ``subpulse``, with ``incident_wavelength``, ``incident_energy``,
        ``time_at_sample``, and ``earliest_arrival``, ``elastic_arrival`` and ``latest_arrival``
        at the detector.

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
    return cli.markdown_table(header, rows)


def elastic_arrival(
    time_at_sample: sc.Variable, incident_wavelength: sc.Variable, L2: sc.Variable
) -> sc.Variable:
    """Return when a neutron that flies through the sample without changing speed arrives.

    Parameters
    ----------
    time_at_sample:
        When each sub-pulse reaches the sample.
    incident_wavelength:
        The wavelength of each sub-pulse, with the dimensions of ``time_at_sample``.
    L2:
        The distance from the sample to the detector.

    Returns
    -------
    :
        The time of arrival at the detector, in the unit of ``time_at_sample``.

    Examples
    --------
    A neutron of 2 Å flies 3 m in 1.5167 ms:

    >>> arrival = elastic_arrival(
    ...     sc.scalar(10.0, unit="ms"), sc.scalar(2.0, unit="Å"), sc.scalar(3.0, unit="m")
    ... )
    >>> print(f"{arrival.value:.4f}")
    11.5167
    """
    flight = tof_kernels.tof_from_wavelength(wavelength=incident_wavelength, Ltotal=L2)
    return time_at_sample + flight.to(unit=time_at_sample.unit)
