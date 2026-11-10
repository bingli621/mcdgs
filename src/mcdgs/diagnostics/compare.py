"""McStas and tof at one monitor, side by side.

:func:`read_comparison` reads one time and wavelength map from both runs, in n/s, and
:func:`plot_comparison` draws it. :func:`read_mcstas_map` reads the map of one McStas run. One
run of each model can feed several monitors: :func:`read_scaled_comparisons` reads them all, tof
scaled at the first and, at the others, with the transmission on top of it.

tof has no flux: each neutron has weight 1. :func:`read_comparison` scales its counts by
``scale``, the n/s each tof neutron stands for. Left out, it is the total of McStas over that of
tof at this monitor, so both have the same total and only their shapes differ. Further down the
beam, pass the scale found where both see the same beam, at the source, for a difference in what
the guides and choppers let through to show. ``transmission`` then scales tof on top of it, for
what the beam loses that tof does not model, such as the guides: tof has none.

:func:`.monitor_bins` gives the bin edges of a monitor's map, read from its ``options``.
"""

from collections.abc import Mapping
from dataclasses import dataclass

import matplotlib.pyplot as plt
import mcstasscript as ms
import numpy as np
import scipp as sc
import tof
from matplotlib.colors import LogNorm, Normalize
from matplotlib.figure import Figure
from mcstasscript.data.data import McStasData
from numpy.typing import NDArray

from mcdgs.tof.reader import read_tof_monitor


@dataclass(frozen=True)
class MonitorComparison:
    """One time and wavelength map in both models, in n/s per bin.

    Each map has the shape (time bins, wavelength bins).

    Parameters
    ----------
    name:
        The name of the monitor, in both models.
    bins:
        The bin edges, ``"toa"`` and ``"wavelength"``, in that order, evenly spaced, in units of
        time and of length.
    mcstas_intensity:
        The intensity of McStas per bin.
    mcstas_error:
        Its error per bin. The bins are independent, so summed over an axis their errors add in
        quadrature.
    tof_intensity:
        The intensity of tof per bin: its counts times ``scale`` and ``transmission``.
    tof_error:
        Its error per bin: the Poisson error of the counts, √N, times ``scale`` and
        ``transmission``.
    ncount:
        Number of McStas neutrons.
    neutrons:
        Number of tof neutrons.
    scale:
        The n/s each tof neutron stands for.
    transmission:
        The fraction of the beam tof is scaled to on top of ``scale``.
    """

    name: str
    bins: Mapping[str, sc.Variable]
    mcstas_intensity: NDArray
    mcstas_error: NDArray
    tof_intensity: NDArray
    tof_error: NDArray
    ncount: float
    neutrons: int
    scale: float
    transmission: float


def read_mcstas_map(
    name: str, mcstas_data: list[McStasData], bins: Mapping[str, sc.Variable]
) -> tuple[NDArray, NDArray]:
    """Read the time and wavelength map of a monitor of one McStas run, in n/s per bin.

    Parameters
    ----------
    name:
        The name of the monitor.
    mcstas_data:
        The monitors of the McStas run, as ``mcstasscript.load_data`` reads them. This one must
        be a 2D map of time against wavelength, ``"time ... lambda ..."``, on ``bins``.
    bins:
        The bin edges of the map, ``"toa"`` and ``"wavelength"``, in that order.

    Returns
    -------
    :
        The intensity and its error per bin, each of shape (time bins, wavelength bins).

    Raises
    ------
    ValueError
        If the map does not have the shape of ``bins``.
    """
    monitor = ms.name_search(name, mcstas_data)
    # McStasScript has the map as (wavelength, time), the y axis first
    intensity = np.asarray(monitor.Intensity, dtype=float).T
    error = np.asarray(monitor.Error, dtype=float).T
    shape = tuple(len(edges) - 1 for edges in bins.values())
    if intensity.shape != shape:
        raise ValueError(
            f"{name}: McStas has a map of {intensity.shape} bins, the bins are {shape}"
        )
    return intensity, error


def read_comparison(
    name: str,
    mcstas_data: list[McStasData],
    tof_result: tof.Result,
    bins: Mapping[str, sc.Variable],
    ncount: float,
    neutrons: int,
    scale: float | None = None,
    transmission: float = 1.0,
) -> MonitorComparison:
    """Read one monitor from a McStas and a tof run, in n/s per bin.

    Parameters
    ----------
    name:
        The name of the monitor, in both models.
    mcstas_data:
        The monitors of the McStas run, as ``mcstasscript.load_data`` reads them. This one must
        be a 2D map of time against wavelength, ``"time ... lambda ..."``, on ``bins``.
    tof_result:
        The tof run.
    bins:
        The bin edges of the McStas map, ``"toa"`` and ``"wavelength"``, in that order.
    ncount:
        Number of McStas neutrons.
    neutrons:
        Number of tof neutrons.
    scale:
        The n/s each tof neutron stands for. If ``None``, the total of McStas over that of tof
        at this monitor.
    transmission:
        A factor on tof on top of ``scale``, for what the beam loses that tof does not model,
        such as the guides. Only with a given ``scale``: with that of this monitor, tof would
        no longer have the total of McStas.

    Returns
    -------
    :
        Both maps, with their errors.

    Raises
    ------
    ValueError
        If the map of McStas does not have the shape of ``bins``, or ``transmission`` is not 1
        without a given ``scale``.
    """
    if scale is None and transmission != 1:
        raise ValueError(f"{name}: a transmission of {transmission} needs a given scale")
    mcstas_intensity, mcstas_error = read_mcstas_map(name, mcstas_data, bins)

    tof_counts = read_tof_monitor(tof_result[name], bins)
    if scale is None:
        scale = float(mcstas_intensity.sum() / tof_counts.sum())

    return MonitorComparison(
        name=name,
        bins=bins,
        mcstas_intensity=mcstas_intensity,
        mcstas_error=mcstas_error,
        tof_intensity=tof_counts * scale * transmission,  # in n/s, times the transmission
        tof_error=np.sqrt(tof_counts) * scale * transmission,  # Poisson, scaled as the counts
        ncount=ncount,
        neutrons=neutrons,
        scale=scale,
        transmission=transmission,
    )


def read_scaled_comparisons(
    mcstas_data: list[McStasData],
    tof_result: tof.Result,
    bins: Mapping[str, Mapping[str, sc.Variable]],
    ncount: float,
    neutrons: int,
    transmission: float | None = None,
) -> dict[str, MonitorComparison]:
    """Read several monitors from a McStas and a tof run, tof scaled at the first one.

    The first monitor is where both models see the same beam, such as the source monitor: there
    tof is scaled to the total of McStas. The others, further down the beam, have that scale and
    ``transmission`` on top of it, as :func:`read_comparison` explains.

    Parameters
    ----------
    mcstas_data:
        The monitors of the McStas run, as ``mcstasscript.load_data`` reads them.
    tof_result:
        The tof run.
    bins:
        The bin edges of each McStas map, ``"toa"`` and ``"wavelength"``, by the name of the
        monitor, in beam order: at least two monitors.
    ncount:
        Number of McStas neutrons.
    neutrons:
        Number of tof neutrons.
    transmission:
        A factor on tof at every monitor after the first, on top of the scale of the first. If
        ``None``, tuned so that tof has the total of McStas at the second monitor.

    Returns
    -------
    :
        Every monitor, by name, in the order of ``bins``.

    Raises
    ------
    ValueError
        If ``bins`` has fewer than two monitors.
    """
    if len(bins) < 2:
        raise ValueError(f"needs a monitor to scale tof at and one after it, got {list(bins)}")
    (first, first_bins), *rest = bins.items()
    at_first = read_comparison(first, mcstas_data, tof_result, first_bins, ncount, neutrons)

    def after_first(
        name: str, edges: Mapping[str, sc.Variable], factor: float
    ) -> MonitorComparison:
        """Read a monitor after the first, with its scale and a transmission."""
        return read_comparison(
            name,
            mcstas_data,
            tof_result,
            edges,
            ncount,
            neutrons,
            scale=at_first.scale,
            transmission=factor,
        )

    if transmission is None:  # tuned: the transmission that gives tof the total of McStas
        second = after_first(*rest[0], 1.0)
        transmission = float(second.mcstas_intensity.sum() / second.tof_intensity.sum())
    return {first: at_first} | {
        name: after_first(name, edges, transmission) for name, edges in rest
    }


def plot_comparison(
    comparison: MonitorComparison,
    title: str,
    *,
    log_decades: int | None = None,
    connect: bool = False,
    log_y: bool = False,
) -> Figure:
    """Draw both models at one monitor, in n/s.

    The top row has the time and wavelength map of each, on one color scale; the bottom row
    the wavelength spectrum and the time profile of both, summed over the other axis, with
    their error bars.

    Parameters
    ----------
    comparison:
        Both maps.
    title:
        The title of the figure.
    log_decades:
        If given, the color scale of the maps is logarithmic, this many decades deep below the
        highest bin, in viridis with the empty bins in white: for a beam whose bins differ by
        orders of magnitude, such as the sub-pulses the choppers cut. Its weakest bins are dark,
        so they stand out on the white. If ``None``, linear from 0, in turbo.
    connect:
        True joins the points of the spectrum and the profile with lines.
    log_y:
        True puts the intensity of the spectrum and the profile on a log axis, where the empty
        bins are left out. False, a linear one from 0.

    Returns
    -------
    :
        The figure.
    """
    c = comparison
    times = c.bins["toa"].to(unit="ms").values
    wavelengths = c.bins["wavelength"].to(unit="Å").values
    time_bin, wavelength_bin = times[1] - times[0], wavelengths[1] - wavelengths[0]
    vmax = max(c.mcstas_intensity.max(), c.tof_intensity.max())  # one color scale for both
    maps = [c.mcstas_intensity.T, c.tof_intensity.T]
    norm: Normalize = Normalize(vmin=0, vmax=vmax)
    cmap = plt.get_cmap("turbo")
    if log_decades is not None:  # the empty bins, masked, in white
        maps = [np.ma.masked_less_equal(m, 0) for m in maps]
        norm = LogNorm(vmin=vmax / 10**log_decades, vmax=vmax, clip=True)
        cmap = plt.get_cmap("viridis").with_extremes(bad="white")

    fig, ((mcstas_ax, tof_ax), (spectrum_ax, profile_ax)) = plt.subplots(
        2, 2, figsize=(12, 9), layout="constrained"
    )
    fig.suptitle(title)

    # the time and wavelength maps
    mcstas_ax.pcolormesh(times, wavelengths, maps[0], cmap=cmap, norm=norm)
    mesh = tof_ax.pcolormesh(times, wavelengths, maps[1], cmap=cmap, norm=norm)
    mcstas_ax.set(title="McStas", xlabel="Time [ms]", ylabel="Wavelength [Å]")
    tof_ax.set(title="tof", xlabel="Time [ms]", ylabel="Wavelength [Å]")
    fig.colorbar(
        mesh,
        ax=[mcstas_ax, tof_ax],
        label=f"Intensity [n/s per {time_bin:.2g} ms and {wavelength_bin:.2g} Å]",
    )

    time_centers = sc.midpoints(c.bins["toa"]).to(unit="ms").values
    wavelength_centers = sc.midpoints(c.bins["wavelength"]).to(unit="Å").values
    errorbar_style = {
        "fmt": ".-" if connect else ".",  # the line only joins the points
        "markersize": 3,
        "linewidth": 0.8,
        "elinewidth": 0.8,
    }
    mcstas_label = f"McStas, ncount = {int(c.ncount):,}"
    tof_label = f"tof, {c.neutrons:,} neutrons × {c.scale:.3g} n/s each"
    if c.transmission != 1:
        tof_label += f"\ntransmission {c.transmission:.3g}"

    # the wavelength spectrum, summed over time, and the time profile, summed over wavelength
    for ax, centers, axis in ((spectrum_ax, wavelength_centers, 0), (profile_ax, time_centers, 1)):
        for intensity, error, label in (
            (c.mcstas_intensity, c.mcstas_error, mcstas_label),
            (c.tof_intensity, c.tof_error, tof_label),
        ):
            ax.errorbar(
                centers,
                intensity.sum(axis=axis),
                yerr=np.sqrt((error**2).sum(axis=axis)),
                label=label,
                **errorbar_style,
            )
        ax.legend()
    spectrum_ax.set(xlabel="Wavelength [Å]", ylabel=f"Intensity [n/s per {wavelength_bin:.2g} Å]")
    profile_ax.set(xlabel="Time [ms]", ylabel=f"Intensity [n/s per {time_bin:.2g} ms]")
    for ax in (spectrum_ax, profile_ax):
        if log_y:  # the empty bins, and the error bars below 0, are left out
            ax.set_yscale("log", nonpositive="mask")
        else:
            ax.set_ylim(bottom=0)
    return fig
