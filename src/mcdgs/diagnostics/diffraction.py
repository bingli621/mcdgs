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

The frame is that of ``sample_pos``, the beam frame of a run from the sample: the origin is the
sample, z along the beam, y up and x to the left of the beam.
"""

import matplotlib.pyplot as plt
import numpy as np
import scipp as sc
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure

from mcdgs.reduction.binning import bin_edges, bins_title
from mcdgs.reduction.inelastic import DEFAULT_ENERGY_TRANSFER_RATIO

Q_STEP = 0.02  # 1/Å, the bin width in Qx and Qz
LOG_DECADES = 4  # the color scale spans this many decades below the highest bin, by default
QY_HALF_WIDTH = 0.5  # 1/Å: the plane keeps the neutrons with |Qy| below this, near horizontal
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
    edges = {}
    for name in ("Qx", "Qz"):
        values = events.coords[name].to(unit="1/Å").values
        start, stop = np.floor(values.min() * 10) / 10, np.ceil(values.max() * 10) / 10
        edges[name] = bin_edges(name, start, stop, step, unit="1/Å")
    plane = events.hist(Qz=edges["Qz"], Qx=edges["Qx"])  # summed over the slice in Qy
    return draw_plane(plane, f"{title}\n{plot_bins(plane, qy_half_width)}", log_decades=log_decades)


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
    plane: sc.DataArray, direct_beam_radius: float = DIRECT_BEAM_RADIUS
) -> float:
    """Return the highest bin of a (Qz, Qx) histogram beyond the direct beam.

    Parameters
    ----------
    plane:
        The histogram, with the bin edges ``Qx`` and ``Qz``, in 1/Å, and any dimensions before
        them.
    direct_beam_radius:
        The |Q| of the direct beam, in 1/Å, left out.

    Returns
    -------
    :
        The highest value of the bins with ``|Q| > direct_beam_radius``.
    """
    qx, qz = (sc.midpoints(plane.coords[name]).to(unit="1/Å").values for name in ("Qx", "Qz"))
    outside = np.hypot(*np.meshgrid(qz, qx, indexing="ij")) > direct_beam_radius
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
