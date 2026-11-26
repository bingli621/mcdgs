"""A rotation scan of a crystal, reduced and drawn, for any instrument.

The rotation scan of each instrument turns the La2Ni7 crystal of :func:`.la2ni7`, with
``IN_PLANE`` in the horizontal plane, by the instrument parameter ``omega`` about the vertical,
and runs the part from ``sample_pos`` to the detector once per angle, from the beam file of its
diffraction diagnostic. This module does the rest:

- :func:`reduce_scan` reduces each angle as :mod:`.diagnostics.diffraction` does, per incident
  neutron and per steradian, turns the momentum transfer into the frame of the crystal with
  :func:`crystal_frame`, and keeps the neutrons near the horizontal plane,
  ``|Qy| < QY_HALF_WIDTH``, in one (Qx, Qz) map per angle.
- :func:`crop` cuts the maps to the bins with neutrons.
- :func:`plot_scan` draws the whole scan, the sum of the maps, and :func:`animate` saves it as a
  GIF that adds one angle per frame, both on the log color scale of :func:`.draw_plane`.
"""

from collections.abc import Callable
from pathlib import Path

import matplotlib.pyplot as plt
import mcstastox
import numpy as np
import scipp as sc
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.artist import Artist
from matplotlib.figure import Figure

from mcdgs.core import rotation
from mcdgs.diagnostics import diffraction as pattern
from mcdgs.mcstas import numeric_parameters
from mcdgs.reduction.binning import bin_edges, bins_title
from mcdgs.reduction.inelastic import normalized_events

IN_PLANE = [2, 1, 0]  # hkl of the reciprocal vector in the horizontal plane at omega = 0
Q_STEP = pattern.Q_STEP
QY_HALF_WIDTH = pattern.QY_HALF_WIDTH
FRAMES_PER_SECOND = 4  # of the GIF, one frame per angle


def crystal_frame(events: sc.DataArray, omega: float) -> sc.DataArray:
    """Turn the momentum transfer of neutrons into the frame of a crystal turned by ``omega``.

    The frame is that of the crystal at ``omega`` = 0, the beam frame of ``sample_pos``, so
    every angle of a scan lands its reflections at the same Q.

    Parameters
    ----------
    events:
        The neutrons of :func:`.events_to_q_and_energy`, with ``Qx``, ``Qy`` and ``Qz``.
    omega:
        The angle of the crystal about the vertical, in deg, as McStas ``ROTATED (0, omega, 0)``.

    Returns
    -------
    :
        The neutrons, with ``Qx``, ``Qy`` and ``Qz`` in the frame of the crystal.

    Examples
    --------
    Turned by 90° about the vertical, the crystal's x is the lab's -z:

    >>> events = sc.DataArray(
    ...     sc.ones(dims=["event"], shape=[1]),
    ...     coords={
    ...         name: sc.array(dims=["event"], values=[value], unit="1/Å")
    ...         for name, value in (("Qx", 0.0), ("Qy", 0.0), ("Qz", -1.0))
    ...     },
    ... )
    >>> turned = crystal_frame(events, 90.0)
    >>> [round(float(turned.coords[n].values[0]), 6) + 0.0 for n in ("Qx", "Qy", "Qz")]
    [1.0, 0.0, 0.0]
    """
    names = ("Qx", "Qy", "Qz")
    unit = events.coords["Qx"].unit
    lab = np.stack([events.coords[name].to(unit=unit).values for name in names])
    turned = rotation(0.0, omega, 0.0) @ lab  # McStas: from the lab frame to the crystal's
    events = events.copy(deep=False)
    for name, values in zip(names, turned):
        events.coords[name] = sc.array(dims=events.dims, values=values, unit=unit)
    return events


def reduce_scan(
    scan: Path,
    pulses: sc.DataArray,
    read_detector: Callable[[Path, int], sc.DataArray],
    pixel_solid_angle: Callable[[sc.DataArray], sc.Variable],
    q_limit: float,
    step: float = Q_STEP,
    qy_half_width: float = QY_HALF_WIDTH,
) -> sc.DataArray:
    """Reduce every angle of a scan into its own (Qx, Qz) map, in the crystal's frame.

    Each angle is read and reduced on its own, so the scan never needs all its neutrons in memory
    at once.

    Parameters
    ----------
    scan:
        The scan folder, with one entry per angle and the angle as the instrument parameter
        ``omega``.
    pulses:
        The sub-pulses of the run to the sample whose beam the scan started from, the same for
        every angle.
    read_detector:
        Reads the neutrons of the detector of an entry of the scan, from 1, by pixel, as the
        ``read_detector`` of the instrument's ``sample_runs``.
    pixel_solid_angle:
        The solid angle of each pixel of what ``read_detector`` reads.
    q_limit:
        The maps span ``-q_limit`` to ``q_limit`` in Qx and Qz, in 1/Å.
    step:
        The bin width in Qx and Qz, in 1/Å.
    qy_half_width:
        Half the width of the slice in Qy, in the crystal's frame, in 1/Å.

    Returns
    -------
    :
        The intensity per bin of :func:`.normalized_events`, per incident neutron and per
        steradian, summed over the slice in Qy, along ``omega``, ``Qz`` and ``Qx``, with the
        angle of each map as the coordinate ``omega``, in deg, and the bin edges ``Qz`` and
        ``Qx``.

    Raises
    ------
    ValueError
        If the scan has no angles.
    """
    edges = {name: bin_edges(name, -q_limit, q_limit, step, unit="1/Å") for name in ("Qz", "Qx")}
    with mcstastox.Read(str(scan)) as data:
        n_entries = data.get_number_of_entries()
    if n_entries == 0:
        raise ValueError(f"{scan}: the scan has no angles")
    maps, omegas = [], []
    for entry in range(1, n_entries + 1):
        with mcstastox.Read(str(scan), entry_number=entry) as data:
            omega = numeric_parameters(data.get_instrument_parameters())["omega"]
        binned = read_detector(scan, entry)
        events = normalized_events(binned, pulses, pixel_solid_angle(binned))
        events = crystal_frame(events, omega)
        qy = events.coords["Qy"].to(unit="1/Å").values
        (dim,) = events.dims
        near = events[dim, np.flatnonzero(np.abs(qy) < qy_half_width)]
        maps.append(near.hist(Qz=edges["Qz"], Qx=edges["Qx"]))
        omegas.append(omega)
        print(f"omega {omega:+6.1f}°: {near.sizes[dim]:9d} neutrons near the plane")
    stack = sc.concat(maps, "omega")
    stack.coords["omega"] = sc.array(dims=["omega"], values=omegas, unit="deg")
    return stack


def crop(plane: sc.DataArray) -> sc.DataArray:
    """Cut a (Qz, Qx) histogram to the bins from the first to the last with neutrons.

    A histogram with dimensions before Qz and Qx, such as ``omega``, is cut to the bins with
    neutrons in any of them.

    Examples
    --------
    >>> plane = sc.DataArray(
    ...     sc.array(dims=["Qz", "Qx"], values=[[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]),
    ...     coords={"Qz": sc.arange("Qz", 4.0), "Qx": sc.arange("Qx", 4.0)},
    ... )
    >>> crop(plane).sizes
    {'Qz': 1, 'Qx': 1}
    """
    values = plane.transpose([*(d for d in plane.dims if d not in ("Qz", "Qx")), "Qz", "Qx"])
    filled = (values.values > 0).reshape(-1, *values.values.shape[-2:]).any(axis=0)
    rows = np.flatnonzero(filled.any(axis=1))
    columns = np.flatnonzero(filled.any(axis=0))
    if rows.size == 0:
        return plane
    return plane["Qz", rows[0] : rows[-1] + 1]["Qx", columns[0] : columns[-1] + 1]


def angles_summed(omegas: np.ndarray, n: int) -> str:
    """Describe the first ``n`` angles of a scan, such as in the title of their sum.

    Examples
    --------
    >>> angles_summed(np.array([-30.0, -28.0, -26.0]), 2)
    'omega from -30.0° to -28.0°, 2 angles summed'
    """
    return f"omega from {omegas[0]:+.1f}° to {omegas[n - 1]:+.1f}°, {n} angles summed"


def plot_scan(
    maps: sc.DataArray,
    title: str,
    qy_half_width: float = QY_HALF_WIDTH,
    log_decades: float = pattern.LOG_DECADES,
) -> Figure:
    """Draw the whole scan: the sum of the maps of every angle, as :func:`.draw_plane`.

    Parameters
    ----------
    maps:
        The maps of :func:`reduce_scan`, along ``omega``, in the order of the scan.
    title:
        The title of the figure; the bins, with the angles summed over, are added below it, as
        :func:`.bins_title`.
    qy_half_width:
        Half the width of the slice in Qy the maps were summed over, in 1/Å.
    log_decades:
        The depth of the color scale, as :func:`.draw_plane`.

    Returns
    -------
    :
        The figure, Qx along ``IN_PLANE`` and Qz along (0 0 1).
    """
    hkl = " ".join(str(i) for i in IN_PLANE)
    summed = maps.sum("omega")
    omegas = maps.coords["omega"].to(unit="deg")
    bins = pattern.plot_bins(summed, qy_half_width)
    return pattern.draw_plane(
        summed,
        f"{title}\n{bins}, {bins_title({'omega': omegas})}",
        xlabel=f"Qx [1/Å], along ({hkl})",
        ylabel="Qz [1/Å], along (0 0 1)",
        log_decades=log_decades,
    )


def animate(
    maps: sc.DataArray,
    path: Path,
    title: str,
    fps: int = FRAMES_PER_SECOND,
    qy_half_width: float = QY_HALF_WIDTH,
    log_decades: float = pattern.LOG_DECADES,
) -> Path:
    """Save the scan as a GIF that adds one angle per frame, with the angles so far in its title.

    Frame ``n`` is the sum of the maps of the first ``n + 1`` angles, so the map fills in as the
    crystal turns, and the last frame is the whole scan, :func:`plot_scan`. Every frame has
    its color scale, so that the frames compare.

    Parameters
    ----------
    maps:
        The maps of :func:`reduce_scan`, along ``omega``, in the order of the scan.
    path:
        The GIF to write.
    title:
        The title of every frame; the bins and the angles so far are added below it.
    fps:
        Frames per second.
    qy_half_width, log_decades:
        As :func:`plot_scan`.

    Returns
    -------
    :
        ``path``.
    """
    omegas = maps.coords["omega"].to(unit="deg").values
    values = maps.transpose(["omega", "Qz", "Qx"]).values
    summed = np.cumsum(values, axis=0)  # frame n: the angles up to n
    # the axes and the color scale of the whole scan
    fig = plot_scan(maps, title, qy_half_width, log_decades)
    (mesh,) = fig.axes[0].collections
    bins = pattern.plot_bins(maps.sum("omega"), qy_half_width)

    def frame(n: int) -> list[Artist]:
        mesh.set_array(np.ma.masked_less_equal(summed[n], 0))
        fig.suptitle(f"{title}\n{bins}, {angles_summed(omegas, n + 1)}")
        return [mesh]

    animation = FuncAnimation(fig, frame, frames=len(omegas))
    animation.save(path, writer=PillowWriter(fps=fps), dpi=100)
    plt.close(fig)
    return path
