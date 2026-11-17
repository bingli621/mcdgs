"""McStas and tof at one monitor, side by side.

Each model gives one time and wavelength map of the monitor, a scipp data array with its bin
edges and its errors: :func:`.read_mcstas_map` that of McStas, in n/s, and :func:`.read_tof_map`
that of tof, in neutrons, both on the bins of the monitor (:func:`.monitor_bins`). The small
pieces work on any such map:

- :func:`mean_coordinate` gives its mean time or wavelength;
- :func:`plot_map` draws it, and :func:`plot_profile` what is left of it summed over an axis,
  such as its wavelength spectrum, with error bars;
- :func:`intensity_label` names its intensity, per bin.

tof has no flux: each neutron has weight 1. Its counts are scaled by ``scale``, the n/s each tof
neutron stands for: the total of McStas over that of tof where both see the same beam, at the
source. There both then have the same total, and only their shapes differ. Further down the
beam, the same scale lets a difference in what the guides and choppers let through show.
``transmission`` then scales tof on top of it, for what the beam loses that tof does not model,
such as the guides: tof has none.

Scaling tof is a few small steps, at each monitor on the same bins: :func:`read_maps` reads the
map of each model, :func:`tof_scale` gives the scale where both see the same beam,
:func:`tuned_transmission` the transmission that gives tof the total of McStas further down, and
:func:`compare_maps` puts both maps side by side in n/s, as a :class:`MonitorComparison`.
:func:`plot_comparison` draws a comparison and :func:`print_means` prints the mean time and
wavelength of each model.
"""

from collections.abc import Mapping
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import scipp as sc
import tof
from matplotlib.axes import Axes
from matplotlib.collections import QuadMesh
from matplotlib.colors import LogNorm, Normalize
from matplotlib.figure import Figure
from mcstasscript.data.data import McStasData
from numpy.typing import NDArray

from mcdgs.mcstas import read_mcstas_map
from mcdgs.tof import read_tof_map

# how a coordinate of a map is shown, by its name: the label of its axis and its unit
AXES: dict[str, tuple[str, str]] = {"toa": ("Time", "ms"), "wavelength": ("Wavelength", "Å")}


def _shown(data: sc.DataArray, dim: str) -> tuple[str, sc.Variable]:
    """Return the axis label and the bin edges of a dimension of a map, as they are shown."""
    name, unit = AXES.get(dim, (dim, str(data.coords[dim].unit)))
    return f"{name} [{unit}]", data.coords[dim].to(unit=unit)


def intensity_label(data: sc.DataArray) -> str:
    """Return the label of the intensity of a map, with the width of its bins.

    Parameters
    ----------
    data:
        A map in n/s per bin, with evenly spaced bin edges as coordinates.

    Returns
    -------
    :
        Such as ``"Intensity [n/s per 0.5 ms and 0.5 Å]"``, with each dimension in the unit it
        is shown in.

    Examples
    --------
    >>> edges = sc.linspace("toa", 0.0, 0.015, 31, unit="s")
    >>> intensity_label(sc.DataArray(sc.zeros(dims=["toa"], shape=[30]), coords={"toa": edges}))
    'Intensity [n/s per 0.5 ms]'
    """
    widths = []
    for dim in data.dims:
        edges = _shown(data, dim)[1]
        widths.append(f"{(edges[1] - edges[0]).value:.2g} {edges.unit}")
    return f"Intensity [n/s per {' and '.join(widths)}]"


def mean_coordinate(data: sc.DataArray, dim: str) -> sc.Variable:
    """Return the mean of a coordinate of a map, weighted by its intensity.

    Parameters
    ----------
    data:
        A map, with bin edges as coordinates.
    dim:
        The dimension, such as ``"toa"`` or ``"wavelength"``.

    Returns
    -------
    :
        The mean of the bin centers of ``dim``, each weighted by the intensity of the map
        summed over its other dimensions, in the unit of the coordinate.

    Examples
    --------
    >>> edges = sc.array(dims=["wavelength"], values=[1.0, 2.0, 3.0], unit="Å")
    >>> counts = sc.array(dims=["wavelength"], values=[1.0, 3.0])
    >>> mean = mean_coordinate(sc.DataArray(counts, coords={"wavelength": edges}), "wavelength")
    >>> float(mean.value), str(mean.unit)
    (2.25, 'Å')
    """
    weights = sc.values(data.data).sum([other for other in data.dims if other != dim])
    return (sc.midpoints(data.coords[dim]) * weights).sum() / weights.sum()


def plot_map(
    ax: Axes, data: sc.DataArray, *, vmax: float | None = None, log_decades: int | None = None
) -> QuadMesh:
    """Draw a two-dimensional map in color, its first dimension along x.

    Parameters
    ----------
    ax:
        The axes to draw on. Their axis labels are set.
    data:
        The map, such as time by wavelength, with bin edges as coordinates.
    vmax:
        The top of the color scale. If ``None``, the highest bin of the map. Give two maps the
        same one to draw them on one color scale.
    log_decades:
        If given, the color scale is logarithmic, this many decades deep below ``vmax``, in
        viridis with the empty bins in white: for a beam whose bins differ by orders of
        magnitude, such as the sub-pulses the choppers cut. Its weakest bins are dark, so they
        stand out on the white. If ``None``, linear from 0, in turbo.

    Returns
    -------
    :
        The mesh, for a color bar.
    """
    (x_label, x_edges), (y_label, y_edges) = (_shown(data, dim) for dim in data.dims)
    values = data.values.T  # pcolormesh has the rows along y
    vmax = float(values.max()) if vmax is None else vmax
    norm: Normalize = Normalize(vmin=0, vmax=vmax)
    cmap = plt.get_cmap("turbo")
    if log_decades is not None:  # the empty bins, masked, in white
        values = np.ma.masked_less_equal(values, 0)
        norm = LogNorm(vmin=vmax / 10**log_decades, vmax=vmax, clip=True)
        cmap = plt.get_cmap("viridis").with_extremes(bad="white")
    mesh = ax.pcolormesh(x_edges.values, y_edges.values, values, cmap=cmap, norm=norm)
    ax.set(xlabel=x_label, ylabel=y_label)
    return mesh


def plot_profile(ax: Axes, data: sc.DataArray, label: str, *, connect: bool = False) -> None:
    """Draw a one-dimensional profile as points with error bars, at the centers of its bins.

    Parameters
    ----------
    ax:
        The axes to draw on. Their axis labels are set.
    data:
        The profile, with bin edges as its coordinate and variances: a map summed over its other
        dimension, such as the wavelength spectrum ``data.sum("toa")`` or the time profile
        ``data.sum("wavelength")``. The errors of the bins summed add in quadrature, as the
        variances do.
    label:
        Its entry in the legend.
    connect:
        True joins the points with a line.
    """
    (dim,) = data.dims
    x_label, edges = _shown(data, dim)
    ax.errorbar(
        sc.midpoints(edges).values,
        data.values,
        yerr=sc.stddevs(data).values,
        label=label,
        fmt=".-" if connect else ".",  # the line only joins the points
        markersize=3,
        linewidth=0.8,
        elinewidth=0.8,
    )
    ax.set(xlabel=x_label, ylabel=intensity_label(data))


@dataclass(frozen=True)
class MonitorComparison:
    """One time and wavelength map in both models, in n/s per bin.

    Parameters
    ----------
    name:
        The name of the monitor, in both models.
    mcstas:
        The map of McStas, as :func:`.read_mcstas_map`: (time bins, wavelength bins), with its bin
        edges, ``"toa"`` and ``"wavelength"``, evenly spaced, and its errors as variances.
    tof:
        The map of tof on the same bins: its counts times ``scale`` and ``transmission``, and
        with them their Poisson errors.
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
    mcstas: sc.DataArray
    tof: sc.DataArray
    ncount: float
    neutrons: int
    scale: float
    transmission: float

    @property
    def bins(self) -> dict[str, sc.Variable]:
        """The bin edges, ``"toa"`` and ``"wavelength"``, in that order."""
        return {dim: self.mcstas.coords[dim] for dim in self.mcstas.dims}

    @property
    def mcstas_intensity(self) -> NDArray:
        """The intensity of McStas per bin, as an array."""
        return self.mcstas.values

    @property
    def mcstas_error(self) -> NDArray:
        """The error of the intensity of McStas per bin, as an array."""
        return sc.stddevs(self.mcstas).values

    @property
    def tof_intensity(self) -> NDArray:
        """The intensity of tof per bin, as an array."""
        return self.tof.values

    @property
    def tof_error(self) -> NDArray:
        """The error of the intensity of tof per bin, as an array."""
        return sc.stddevs(self.tof).values


def read_maps(
    name: str,
    mcstas_data: list[McStasData],
    tof_result: tof.Result,
    bins: Mapping[str, sc.Variable],
) -> tuple[sc.DataArray, sc.DataArray]:
    """Read one monitor of a McStas and a tof run on the same bins.

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

    Returns
    -------
    :
        The map of McStas in n/s, as :func:`.read_mcstas_map`, and that of tof in neutrons, as
        :func:`.read_tof_map`, both with their errors.

    Raises
    ------
    ValueError
        If the map of McStas does not have the shape of ``bins``.
    """
    return read_mcstas_map(name, mcstas_data, bins), read_tof_map(name, tof_result, bins)


def tof_scale(mcstas_map: sc.DataArray, tof_counts: sc.DataArray) -> float:
    """Return the n/s each tof neutron stands for: the total of McStas over that of tof.

    Set it where both models see the same beam, such as the source monitor, and keep it further
    down the beam, where a new one would hide what the guides and choppers let through.

    Parameters
    ----------
    mcstas_map:
        The map of McStas, in n/s.
    tof_counts:
        The map of tof on the same bins, in neutrons.

    Returns
    -------
    :
        The scale, in n/s per tof neutron.

    Examples
    --------
    >>> mcstas_map = sc.DataArray(sc.array(dims=["toa"], values=[2.0, 4.0], unit="1/s"))
    >>> tof_counts = sc.DataArray(sc.array(dims=["toa"], values=[1.0, 2.0]))
    >>> tof_scale(mcstas_map, tof_counts)
    2.0
    """
    return float(mcstas_map.sum().value / tof_counts.sum().value)


def tuned_transmission(mcstas_map: sc.DataArray, tof_counts: sc.DataArray, scale: float) -> float:
    """Return the factor on top of ``scale`` that gives tof the total of McStas at a monitor.

    Parameters
    ----------
    mcstas_map:
        The map of McStas, in n/s.
    tof_counts:
        The map of tof on the same bins, in neutrons.
    scale:
        The n/s each tof neutron stands for, from :func:`tof_scale` at an earlier monitor.

    Returns
    -------
    :
        The transmission: the total of McStas over that of tof times ``scale``.

    Examples
    --------
    >>> mcstas_map = sc.DataArray(sc.array(dims=["toa"], values=[1.0, 2.0], unit="1/s"))
    >>> tof_counts = sc.DataArray(sc.array(dims=["toa"], values=[1.0, 2.0]))
    >>> tuned_transmission(mcstas_map, tof_counts, scale=4.0)
    0.25
    """
    return tof_scale(mcstas_map, tof_counts) / scale


def compare_maps(
    name: str,
    mcstas_map: sc.DataArray,
    tof_counts: sc.DataArray,
    ncount: float,
    neutrons: int,
    scale: float,
    transmission: float = 1.0,
) -> MonitorComparison:
    """Put both maps of one monitor side by side, tof scaled to n/s.

    Parameters
    ----------
    name:
        The name of the monitor, in both models.
    mcstas_map:
        The map of McStas, in n/s, as :func:`read_maps` reads it.
    tof_counts:
        The map of tof on the same bins, in neutrons.
    ncount:
        Number of McStas neutrons.
    neutrons:
        Number of tof neutrons.
    scale:
        The n/s each tof neutron stands for, from :func:`tof_scale`.
    transmission:
        A factor on tof on top of ``scale``, for what the beam loses that tof does not model,
        such as the guides; from :func:`tuned_transmission` to give tof the total of McStas.

    Returns
    -------
    :
        Both maps, tof times ``scale × transmission``, its Poisson errors scaled along.
    """
    tof_map = tof_counts * (scale * transmission)
    tof_map.unit = "1/s"
    return MonitorComparison(
        name=name,
        mcstas=mcstas_map,
        tof=tof_map,
        ncount=ncount,
        neutrons=neutrons,
        scale=scale,
        transmission=transmission,
    )


def print_means(comparison: MonitorComparison) -> None:
    """Print the mean time and wavelength of each model, weighted by its intensity.

    Parameters
    ----------
    comparison:
        Both maps.
    """
    for label, data in (("McStas", comparison.mcstas), ("tof", comparison.tof)):
        mean_time = mean_coordinate(data, "toa").to(unit="ms").value
        mean_wavelength = mean_coordinate(data, "wavelength").to(unit="Å").value
        print(f"{label:>6}: mean time {mean_time:.3f} ms, mean wavelength {mean_wavelength:.3f} Å")


def plot_comparison(
    comparison: MonitorComparison,
    title: str,
    *,
    log_decades: int | None = None,
    connect: bool = False,
    log_y: bool = False,
) -> Figure:
    """Draw both models at one monitor, in n/s.

    The top row has the time and wavelength map of each, on one color scale (:func:`plot_map`);
    the bottom row the wavelength spectrum and the time profile of both, summed over the other
    axis, with their error bars (:func:`plot_profile`).

    Parameters
    ----------
    comparison:
        Both maps.
    title:
        The title of the figure.
    log_decades:
        If given, the color scale of the maps is logarithmic, this many decades deep below the
        highest bin, as in :func:`plot_map`. If ``None``, linear from 0.
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
    fig, ((mcstas_ax, tof_ax), (spectrum_ax, profile_ax)) = plt.subplots(
        2, 2, figsize=(12, 9), layout="constrained"
    )
    fig.suptitle(title)

    # the time and wavelength maps, on one color scale
    vmax = float(max(c.mcstas.values.max(), c.tof.values.max()))
    plot_map(mcstas_ax, c.mcstas, vmax=vmax, log_decades=log_decades)
    mesh = plot_map(tof_ax, c.tof, vmax=vmax, log_decades=log_decades)
    mcstas_ax.set_title("McStas")
    tof_ax.set_title("tof")
    fig.colorbar(mesh, ax=[mcstas_ax, tof_ax], label=intensity_label(c.mcstas))

    mcstas_label = f"McStas, ncount = {int(c.ncount):,}"
    tof_label = f"tof, {c.neutrons:,} neutrons × {c.scale:.3g} n/s each"
    if c.transmission != 1:
        tof_label += f"\ntransmission {c.transmission:.3g}"

    # the wavelength spectrum, summed over time, and the time profile, summed over wavelength
    time, wavelength = c.mcstas.dims
    for ax, summed in ((spectrum_ax, time), (profile_ax, wavelength)):
        plot_profile(ax, c.mcstas.sum(summed), mcstas_label, connect=connect)
        plot_profile(ax, c.tof.sum(summed), tof_label, connect=connect)
        ax.legend()
        if log_y:  # the empty bins, and the error bars below 0, are left out
            ax.set_yscale("log", nonpositive="mask")
        else:
            ax.set_ylim(bottom=0)
    return fig
