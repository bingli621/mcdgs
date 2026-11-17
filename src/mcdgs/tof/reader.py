"""Read the neutrons of a tof run.

:func:`read_tof_map` counts the neutrons at a component of a tof run, a monitor or a chopper,
on bins of any of their coordinates, such as the time of arrival and the wavelength.
:func:`.monitor_bins` gives those bins from the axes of a ``Monitor_nD``, so that tof is read on
the bins of McStas, as :func:`.read_mcstas_map` reads McStas.
"""

from collections.abc import Mapping

import scipp as sc
import tof


def read_tof_map(
    name: str,
    tof_result: tof.Result,
    bins: Mapping[str, sc.Variable],
    *,
    transmitted: bool = True,
) -> sc.DataArray:
    """Count the neutrons of the first pulse at a component of a tof run, on any bins.

    Parameters
    ----------
    name:
        The name of the component: a monitor, a ``tof.Detector``, or a chopper, a
        ``tof.Chopper``, the only components a tof run has a reading of.
    tof_result:
        The tof run.
    bins:
        The bin edges of each axis, by the name of a coordinate of the neutrons, such as
        ``"toa"`` or ``"wavelength"``, in any unit that converts to the coordinate's.
    transmitted:
        True counts the neutrons the component lets through: of a chopper, those that pass its
        slits, as a McStas monitor just after it would see them. False counts all that reach
        it, also those a chopper blocks. A monitor blocks none: both count the same.

    Returns
    -------
    :
        The neutrons per bin, in counts, with their Poisson variances: N for a bin of N
        neutrons, an error of √N. It has one dimension per entry of ``bins``, in their order,
        and ``bins`` as its coordinates, in the units given.
    """
    neutrons = tof_result[name].data["pulse", 0]
    blocked = neutrons.masks["blocked_by_others"]
    if transmitted and "blocked_by_me" in neutrons.masks:  # a chopper
        blocked = blocked | neutrons.masks["blocked_by_me"]
    # hist leaves out the masked neutrons itself: the masks go, once applied as wanted
    counted = neutrons[~blocked].drop_masks(list(neutrons.masks))
    # scipp bins only on edges in the unit of the coordinate
    edges = {dim: edge.to(unit=counted.coords[dim].unit) for dim, edge in bins.items()}
    counts = counted.hist(edges).values
    return sc.DataArray(
        sc.array(dims=list(bins), values=counts, variances=counts, unit="counts"),
        coords=dict(bins),
    )
