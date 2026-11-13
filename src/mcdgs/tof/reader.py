"""Read the neutrons of a tof run.

:func:`read_tof_monitor` counts the neutrons that reach a monitor on bins of any of their
coordinates, such as the time of arrival and the wavelength. :func:`monitor_bins` gives those
bins from a ``Monitor_nD``, so that tof is read on the bins of McStas.
"""

from collections.abc import Mapping

import scipp as sc
import tof
from numpy.typing import NDArray

from mcdgs.core.components import Monitor_nD

# the Monitor_nD variables tof has, by McStas name: the tof coordinate and the unit of McStas
TOF_COORDINATES: dict[str, tuple[str, str]] = {"time": ("toa", "s"), "lambda": ("wavelength", "Å")}


def monitor_bins(monitor: Monitor_nD) -> dict[str, sc.Variable]:
    """Return the bin edges of each axis of a monitor that tof has, as tof's coordinate.

    Parameters
    ----------
    monitor:
        A monitor whose ``options`` give its axes, see :attr:`.Monitor_nD.bin_limits`.

    Returns
    -------
    :
        The edges by tof coordinate, such as ``"toa"`` in s and ``"wavelength"`` in Å, in the
        order of ``options``: that of the axes of the map of McStas. An axis tof does not have,
        one not in ``TOF_COORDINATES``, is left out.

    Examples
    --------
    >>> bins = monitor_bins(Monitor_nD(options="lambda limits=[0.5 5] bins=50"))
    >>> list(bins), str(bins["wavelength"].unit)
    (['wavelength'], 'Å')
    """
    bins = {}
    for variable, (start, stop, n_bins) in monitor.bin_limits.items():
        if variable in TOF_COORDINATES:
            coordinate, unit = TOF_COORDINATES[variable]
            bins[coordinate] = sc.linspace(coordinate, start, stop, n_bins + 1, unit=unit)
    return bins


def read_tof_monitor(monitor: tof.DetectorReading, bins: Mapping[str, sc.Variable]) -> NDArray:  # type: ignore
    """Count the neutrons of the first pulse that reach a tof monitor, on any bins.

    Parameters
    ----------
    monitor:
        The reading of a monitor, such as ``result["sample_monitor"]`` of a tof run.
    bins:
        The bin edges of each axis, by the name of a coordinate of the neutrons, such as
        ``"toa"`` or ``"wavelength"``, in any unit that converts to the coordinate's.

    Returns
    -------
    :
        The neutrons per bin, with one axis per entry of ``bins``, in their order.
    """
    neutrons = monitor.data["pulse", 0]
    reached = neutrons[~neutrons.masks["blocked_by_others"]]
    # scipp bins only on edges in the unit of the coordinate
    edges = {name: edge.to(unit=reached.coords[name].unit) for name, edge in bins.items()}
    return reached.hist(edges).values
