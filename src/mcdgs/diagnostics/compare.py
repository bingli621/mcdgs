"""McStas and tof at one monitor, side by side.

:func:`read_comparison` reads one time and wavelength map from both runs, in n/s, and
:func:`plot_comparison` draws it. One run of each model can feed several monitors.

tof has no flux: each neutron has weight 1. :func:`read_comparison` scales its counts by
``scale``, the n/s each tof neutron stands for. Left out, it is the total of McStas over that of
tof at this monitor, so both have the same total and only their shapes differ. Further down the
beam, pass the scale found where both see the same beam, at the source, for a difference in what
the guides and choppers let through to show.
"""

from collections.abc import Mapping
from dataclasses import dataclass

import matplotlib.pyplot as plt
import mcstasscript as ms
import numpy as np
import scipp as sc
import tof
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
        The intensity of tof per bin: its counts times ``scale``.
    tof_error:
        Its error per bin: the Poisson error of the counts, √N, times ``scale``.
    ncount:
        Number of McStas neutrons.
    neutrons:
        Number of tof neutrons.
    scale:
        The n/s each tof neutron stands for.
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


def read_comparison(
    name: str,
    mcstas_data: list[McStasData],
    tof_result: tof.Result,
    bins: Mapping[str, sc.Variable],
    ncount: float,
    neutrons: int,
    scale: float | None = None,
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

    Returns
    -------
    :
        Both maps, with their errors.

    Raises
    ------
    ValueError
        If the map of McStas does not have the shape of ``bins``.
    """
    monitor = ms.name_search(name, mcstas_data)
    # McStasScript has the map as (wavelength, time), the y axis first
    mcstas_intensity = np.asarray(monitor.Intensity, dtype=float).T
    mcstas_error = np.asarray(monitor.Error, dtype=float).T
    shape = tuple(len(edges) - 1 for edges in bins.values())
    if mcstas_intensity.shape != shape:
        raise ValueError(
            f"{name}: McStas has a map of {mcstas_intensity.shape} bins, the bins are {shape}"
        )

    tof_counts = read_tof_monitor(tof_result[name], bins)
    if scale is None:
        scale = float(mcstas_intensity.sum() / tof_counts.sum())

    return MonitorComparison(
        name=name,
        bins=bins,
        mcstas_intensity=mcstas_intensity,
        mcstas_error=mcstas_error,
        tof_intensity=tof_counts * scale,
        tof_error=np.sqrt(tof_counts) * scale,  # Poisson, scaled as the counts
        ncount=ncount,
        neutrons=neutrons,
        scale=scale,
    )


def plot_comparison(comparison: MonitorComparison, title: str) -> Figure:
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

    fig, ((mcstas_ax, tof_ax), (spectrum_ax, profile_ax)) = plt.subplots(
        2, 2, figsize=(12, 9), layout="constrained"
    )
    fig.suptitle(title)

    # the time and wavelength maps
    mcstas_ax.pcolormesh(times, wavelengths, c.mcstas_intensity.T, cmap="turbo", vmin=0, vmax=vmax)
    mesh = tof_ax.pcolormesh(times, wavelengths, c.tof_intensity.T, cmap="turbo", vmin=0, vmax=vmax)
    mcstas_ax.set(title="McStas", xlabel="Time [ms]", ylabel="Wavelength [Å]")
    tof_ax.set(title="tof", xlabel="Time [ms]", ylabel="Wavelength [Å]")
    fig.colorbar(
        mesh,
        ax=[mcstas_ax, tof_ax],
        label=f"Intensity [n/s per {time_bin:.2g} ms and {wavelength_bin:.2g} Å]",
    )

    time_centers = sc.midpoints(c.bins["toa"]).to(unit="ms").values
    wavelength_centers = sc.midpoints(c.bins["wavelength"]).to(unit="Å").values
    errorbar_style = {"fmt": ".", "markersize": 3, "elinewidth": 0.8}
    mcstas_label = f"McStas, ncount = {int(c.ncount):,}"
    tof_label = f"tof, {c.neutrons:,} neutrons × {c.scale:.3g} n/s each"

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
    spectrum_ax.set(
        xlabel="Wavelength [Å]", ylabel=f"Intensity [n/s per {wavelength_bin:.2g} Å]", ylim=0
    )
    profile_ax.set(xlabel="Time [ms]", ylabel=f"Intensity [n/s per {time_bin:.2g} ms]", ylim=0)
    return fig
