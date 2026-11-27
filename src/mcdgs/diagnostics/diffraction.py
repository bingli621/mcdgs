"""Diffraction from a crystal, reduced and drawn, for any instrument.

The diffraction diagnostic of each instrument runs McStas, finds its sub-pulses, reads its
detector and converts every neutron in the arrival window of its sub-pulse to its momentum
transfer (Qx, Qy, Qz) with :func:`.normalized_events`, whatever its energy transfer: the
inelastic scattering is weak next to the Bragg peaks, so the plane is a diffraction pattern.
This module does the rest:

- :func:`print_subpulse_table` prints the incident energy, ki and time at the sample of each
  sub-pulse, and the fraction of its neutrons that reach the detector.
- :func:`plot_plane` draws the intensity in (Qx, Qz) of every sub-pulse, summed over a slice of
  Qy near the horizontal plane, ``|Qy| < QY_HALF_WIDTH``, with :func:`draw_plane`, on a log
  color scale, and the bins of :func:`plot_bins` under the title. Each sub-pulse has its own
  incident wavevector ki, so its elastic neutrons lie near its own circle, of radius ki through
  the origin, where the detector covers it.
- :func:`plot_vertical_planes` draws the intensity in (Qy, Qz) of each sub-pulse on its own
  panel, summed over ``|Qx| < QX_HALF_WIDTH``, on one log color scale for all.

The frame is that of ``sample_pos``, the beam frame of a run from the sample: the origin is the
sample, z along the beam, y up and x to the left of the beam.
"""

import matplotlib.pyplot as plt
import numpy as np
import scipp as sc
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure

from mcdgs.reduction.binning import bin_edges, bins_title
from mcdgs.reduction.inelastic import DEFAULT_ENERGY_TRANSFER_RATIO

Q_STEP = 0.02  # 1/Å, the bin width in Qx and Qz
LOG_DECADES = 4  # the color scale spans this many decades below the highest bin, by default
QY_HALF_WIDTH = 0.5  # 1/Å: the plane keeps the neutrons with |Qy| below this, near horizontal
QX_HALF_WIDTH = 5.0  # 1/Å: the (Qy, Qz) maps sum the neutrons with |Qx| below this
VMAX_FRACTION = 0.1  # the (Qy, Qz) maps put the top of their scale at this part of the highest bin
PANEL_COLUMNS = 3  # the (Qy, Qz) maps of the sub-pulses are drawn this many to a row
DIRECT_BEAM_RADIUS = 0.2  # 1/Å: the direct beam, far above the rest, is left out of the scale


def print_subpulse_table(pulses: sc.DataArray, events: sc.DataArray) -> None:
    """Print the incident energy, ki and time at the sample of each sub-pulse.

    Parameters
    ----------
    pulses:
        The sub-pulses, with their flux, ``time_at_sample``, ``incident_wavelength`` and
        ``incident_energy``.
    events:
        The neutrons of :func:`.normalized_events`.
    """
    print(
        "  sub-pulse  t_sample [ms]  λi [Å]  Ei [meV]  ki [1/Å]  flux [n/s]  neutrons"
        "  on the detector"
    )
    subpulse = events.coords["subpulse"].values
    n = pulses.sizes["subpulse"]
    counts = np.bincount(subpulse, minlength=n)
    # the fraction of the incident neutrons that reach the detector, in a window
    fractions = np.bincount(
        subpulse, events.values * events.coords["solid_angle"].values, minlength=n
    )
    for i, (t, wavelength, ei, flux, count, fraction) in enumerate(
        zip(
            pulses.coords["time_at_sample"].to(unit="ms").values,
            pulses.coords["incident_wavelength"].to(unit="Å").values,
            pulses.coords["incident_energy"].to(unit="meV").values,
            pulses.values,
            counts,
            fractions,
        )
    ):
        k = 2 * np.pi / wavelength
        print(
            f"  {i:9d}  {t:13.3f}  {wavelength:6.4f}  {ei:8.3f}  {k:8.4f}  {flux:10.4g}"
            f"  {count:8d}  {fraction:15.4g}"
        )


def plot_plane(
    events: sc.DataArray,
    title: str,
    step: float = Q_STEP,
    qy_half_width: float = QY_HALF_WIDTH,
    log_decades: float = LOG_DECADES,
) -> Figure:
    """Draw the intensity in (Qx, Qz) of every sub-pulse, near the horizontal plane.

    Only the neutrons with ``|Qy| < qy_half_width`` count, summed over Qy: the detector also
    sees reflections well above and below the plane, which a sum over all Qy would project
    onto it.

    It is drawn on the log color scale of :func:`draw_plane`, with the bins of
    :func:`plot_bins` under the title.

    Parameters
    ----------
    events:
        The neutrons of :func:`.events_to_q_and_energy`.
    title:
        The title of the figure.
    step:
        The bin width in Qx and Qz, in 1/Å.
    qy_half_width:
        Half the width of the slice in Qy, in 1/Å.
    log_decades:
        How many decades the color scale spans below its top, as :func:`draw_plane`.

    Returns
    -------
    :
        The figure.
    """
    qy = events.coords["Qy"].to(unit="1/Å").values
    (dim,) = events.dims
    events = events[dim, np.flatnonzero(np.abs(qy) < qy_half_width)]
    edges = {name: covering_edges(events, name, step) for name in ("Qx", "Qz")}
    plane = events.hist(Qz=edges["Qz"], Qx=edges["Qx"])  # summed over the slice in Qy
    return draw_plane(plane, f"{title}\n{plot_bins(plane, qy_half_width)}", log_decades=log_decades)


def covering_edges(events: sc.DataArray, name: str, step: float = Q_STEP) -> sc.Variable:
    """Return bin edges in 1/Å that cover every neutron in one component of Q.

    The edges start and stop on the tenths of 1/Å just outside the neutrons.

    Parameters
    ----------
    events:
        The neutrons, with the coordinate ``name``.
    name:
        The component of Q, such as ``"Qx"``.
    step:
        The bin width, in 1/Å.

    Returns
    -------
    :
        The edges, along ``name``.
    """
    values = events.coords[name].to(unit="1/Å").values
    start, stop = np.floor(values.min() * 10) / 10, np.ceil(values.max() * 10) / 10
    return bin_edges(name, start, stop, step, unit="1/Å")


def plot_vertical_planes(
    events: sc.DataArray,
    pulses: sc.DataArray,
    title: str,
    step: float = Q_STEP,
    qx_half_width: float = QX_HALF_WIDTH,
    log_decades: float = LOG_DECADES,
    vmax_fraction: float = VMAX_FRACTION,
    columns: int = PANEL_COLUMNS,
) -> Figure:
    """Draw the intensity in (Qy, Qz) of each sub-pulse on its own panel.

    Only the neutrons with ``|Qx| < qx_half_width`` count, summed over Qx. All the panels share
    the axes and one log color scale, in turbo, with the empty bins in white. Its bottom is
    ``log_decades`` below the highest bin of any sub-pulse beyond the direct beam, as
    :func:`draw_plane`; its top is ``vmax_fraction`` of that bin, so that the Bragg peaks, not
    only the one bin, reach the red end, and the bins above it take the color of the top. The
    bins are written under the title, as :func:`.bins_title`. Qz runs across, along the
    beam, and Qy up.

    Parameters
    ----------
    events:
        The neutrons of :func:`.normalized_events`, with the index ``subpulse`` of the
        sub-pulse of each.
    pulses:
        The sub-pulses, with their ``incident_wavelength`` and ``incident_energy``, to label the
        panels.
    title:
        The title of the figure.
    step:
        The bin width in Qy and Qz, in 1/Å.
    qx_half_width:
        Half the width of the slice in Qx, in 1/Å.
    log_decades:
        How many decades below the highest bin the color scale starts.
    vmax_fraction:
        The top of the color scale, as a part of the highest bin.
    columns:
        How many panels to a row.

    Returns
    -------
    :
        The figure.
    """
    qx = events.coords["Qx"].to(unit="1/Å").values
    (dim,) = events.dims
    events = events[dim, np.flatnonzero(np.abs(qx) < qx_half_width)]
    edges = {name: covering_edges(events, name, step) for name in ("Qy", "Qz")}
    n = pulses.sizes["subpulse"]
    # every sub-pulse, also one without a neutron in the slice, summed over the slice in Qx
    index = events.coords["subpulse"]
    groups = sc.arange("subpulse", n, unit=index.unit, dtype=index.dtype)
    planes = events.group(groups).hist(Qy=edges["Qy"], Qz=edges["Qz"])
    highest = beyond_direct_beam_max(planes, axes=("Qy", "Qz"))
    norm = LogNorm(vmin=highest / 10**log_decades, vmax=highest * vmax_fraction, clip=True)
    cmap = plt.get_cmap("turbo").with_extremes(bad="white")

    rows = -(-n // columns)  # ceiling division
    # panels 5 inches wide, as high as the aspect of the map, plus room for their titles
    span = {name: float(edges[name].values[-1] - edges[name].values[0]) for name in edges}
    row_height = 5.0 * span["Qy"] / span["Qz"] + 0.5
    fig, axes = plt.subplots(
        rows,
        columns,
        figsize=(6.0 * columns, row_height * rows + 1.2),
        sharex=True,
        sharey=True,
        squeeze=False,
        layout="constrained",
    )
    bins = bins_title(
        {
            "Qz": edges["Qz"],
            "Qy": edges["Qy"],
            "Qx": bin_edges("Qx", -qx_half_width, qx_half_width, unit="1/Å"),
            "E/Ei": bin_edges("E/Ei", *DEFAULT_ENERGY_TRANSFER_RATIO),
        }
    )
    fig.suptitle(f"{title}\n{bins}")
    wavelengths = pulses.coords["incident_wavelength"].to(unit="Å").values
    energies = pulses.coords["incident_energy"].to(unit="meV").values
    for i, ax in enumerate(axes.flat):
        if i >= n:
            ax.set_visible(False)  # the panels left over in the last row
            continue
        ax.pcolormesh(
            edges["Qz"].values,
            edges["Qy"].values,
            np.ma.masked_less_equal(planes["subpulse", i].values, 0),  # empty bins in white
            cmap=cmap,
            norm=norm,
        )
        ax.set(aspect="equal")
        ax.set_title(
            f"sub-pulse {i}\nλi = {wavelengths[i]:.3f} Å, Ei = {energies[i]:.2f} meV",
            fontsize="small",
        )
        ax.grid(True, alpha=0.3)
    for ax in axes[-1]:
        ax.set_xlabel("Qz [1/Å], along the beam")
    for ax in axes[:, 0]:
        ax.set_ylabel("Qy [1/Å], up")
    fig.colorbar(
        ScalarMappable(norm=norm, cmap=cmap),  # the scale that all the panels share
        ax=axes,
        shrink=0.6,
        label=f"Intensity [{planes.unit} per ({step:g} 1/Å)²]",
        extend="max",
    )
    return fig


def plot_bins(
    plane: sc.DataArray,
    qy_half_width: float = QY_HALF_WIDTH,
    window: tuple[float, float] = DEFAULT_ENERGY_TRANSFER_RATIO,
) -> str:
    """Write the bins of a (Qz, Qx) plane for its title, with :func:`.bins_title`.

    Qx and Qz are binned; the plane is summed over the slice in Qy and over the arrival window of
    each sub-pulse in E / Ei.

    Parameters
    ----------
    plane:
        The histogram, with the bin edges ``Qx`` and ``Qz``.
    qy_half_width:
        Half the width of the slice in Qy, in 1/Å.
    window:
        The smallest and the largest ``E / Ei`` of the arrival windows, as
        :func:`.arrival_window`.

    Returns
    -------
    :
        The bins, as ``Qx=(start, stop, step) [1/Å]``.
    """
    return bins_title(
        {
            "Qx": plane.coords["Qx"],
            "Qz": plane.coords["Qz"],
            "Qy": bin_edges("Qy", -qy_half_width, qy_half_width, unit="1/Å"),
            "E/Ei": bin_edges("E/Ei", *window),
        }
    )


def beyond_direct_beam_max(
    plane: sc.DataArray,
    direct_beam_radius: float = DIRECT_BEAM_RADIUS,
    axes: tuple[str, str] = ("Qz", "Qx"),
) -> float:
    """Return the highest bin of a histogram in two components of Q beyond the direct beam.

    Parameters
    ----------
    plane:
        The histogram, with the bin edges of ``axes``, in 1/Å, as its last two dimensions in
        that order, and any dimensions before them.
    direct_beam_radius:
        The |Q| in the plane of ``axes`` of the direct beam, in 1/Å, left out.
    axes:
        The two components of Q of the histogram, such as ``("Qz", "Qx")``.

    Returns
    -------
    :
        The highest value of the bins with ``|Q| > direct_beam_radius``.
    """
    first, second = (sc.midpoints(plane.coords[name]).to(unit="1/Å").values for name in axes)
    outside = np.hypot(*np.meshgrid(first, second, indexing="ij")) > direct_beam_radius
    return float(plane.values[..., outside].max())


def draw_plane(
    plane: sc.DataArray,
    title: str,
    xlabel: str = "Qx [1/Å]",
    ylabel: str = "Qz [1/Å], along the beam",
    direct_beam_radius: float = DIRECT_BEAM_RADIUS,
    vmax: float | None = None,
    log_decades: float = LOG_DECADES,
) -> Figure:
    """Draw a histogram in (Qx, Qz) on a log color scale.

    The scale is turbo, ``log_decades`` deep below the highest bin beyond the direct beam,
    ``|Q| > direct_beam_radius``, with the empty bins in white: the direct beam that passes the
    beamstop, a few bins far above all others, would leave the rest at the bottom of the scale,
    so it takes the color of the top instead. The axes have the same scale, so that circles look
    round.

    Parameters
    ----------
    plane:
        The histogram, with the bin edges ``Qx`` and ``Qz``, in 1/Å, evenly spaced.
    title:
        The title of the figure.
    xlabel, ylabel:
        The labels of the axes.
    direct_beam_radius:
        The |Q| of the direct beam, in 1/Å, that the color scale leaves out.
    vmax:
        The top of the color scale. If ``None``, the highest bin beyond the direct beam.
    log_decades:
        How many decades the color scale spans below its top. Fewer put the weak bins lower
        on the scale, darker.

    Returns
    -------
    :
        The figure.
    """
    step = float(np.diff(plane.coords["Qx"].values)[0])
    if vmax is None:
        vmax = beyond_direct_beam_max(plane, direct_beam_radius)

    fig, ax = plt.subplots(figsize=(9, 8), layout="constrained")
    fig.suptitle(title)
    mesh = ax.pcolormesh(
        plane.coords["Qx"].values,
        plane.coords["Qz"].values,
        np.ma.masked_less_equal(plane.values, 0),  # empty bins in white
        cmap=plt.get_cmap("turbo").with_extremes(bad="white"),
        norm=LogNorm(vmin=vmax / 10**log_decades, vmax=vmax, clip=True),
    )
    ax.set(xlabel=xlabel, ylabel=ylabel, aspect="equal")
    ax.grid(True, alpha=0.3)
    fig.colorbar(
        mesh,
        ax=ax,
        shrink=0.8,
        label=f"Intensity [{plane.unit} per ({step:g} 1/Å)²]",
        extend="max",
    )
    return fig
