"""Read the neutrons of a tof run.

:func:`read_tof_monitor` counts the neutrons that reach a monitor on bins of any of their
coordinates, such as the time of arrival and the wavelength.
"""

from collections.abc import Mapping

import scipp as sc
import tof
from numpy.typing import NDArray


def read_tof_monitor(monitor: tof.DetectorReading, bins: Mapping[str, sc.Variable]) -> NDArray:
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
