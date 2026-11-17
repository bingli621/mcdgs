"""Read the neutrons of a tof run.

:func:`read_tof_monitor` counts the neutrons that reach a monitor on bins of any of their
coordinates, such as the time of arrival and the wavelength. :func:`.monitor_bins` gives those
bins from the axes of a ``Monitor_nD``, so that tof is read on the bins of McStas.
"""

from collections.abc import Mapping

import scipp as sc
import tof


def read_tof_monitor(  # type: ignore
    monitor: tof.DetectorReading, bins: Mapping[str, sc.Variable]
) -> sc.DataArray:
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
        The neutrons per bin, in counts, with their Poisson variances: N for a bin of N
        neutrons, an error of √N. It has one dimension per entry of ``bins``, in their order,
        and ``bins`` as its coordinates, in the units given.
    """
    neutrons = monitor.data["pulse", 0]
    reached = neutrons[~neutrons.masks["blocked_by_others"]]
    # scipp bins only on edges in the unit of the coordinate
    edges = {name: edge.to(unit=reached.coords[name].unit) for name, edge in bins.items()}
    counts = reached.hist(edges).values
    return sc.DataArray(
        sc.array(dims=list(bins), values=counts, variances=counts, unit="counts"),
        coords=dict(bins),
    )
