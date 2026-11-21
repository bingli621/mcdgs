"""Read the neutrons of a tof run.

:func:`read_tof_maps` counts the neutrons at components of a tof run, monitors or choppers, on
bins of any of their coordinates, such as the time of arrival and the wavelength.
:func:`.monitor_bins` gives those bins from the axes of a ``Monitor_nD``, so that tof is read on
the bins of McStas, as :func:`.read_mcstas_maps` reads McStas.
"""

from collections.abc import Mapping

import scipp as sc
import tof


def read_tof_maps(
    tof_run: tof.Result,
    bins: Mapping[str, Mapping[str, sc.Variable]],
    *,
    transmitted: bool = True,
) -> dict[str, sc.DataArray]:
    """Count the neutrons of the first pulse at components of a tof run, on any bins.

    Parameters
    ----------
    tof_run:
        The tof run, what ``model.run()`` returns.
    bins:
        The bin edges of each component, by its name, and of each axis, by the name of a
        coordinate of the neutrons, such as ``"toa"`` or ``"wavelength"``, in any unit that
        converts to the coordinate's. A component is a monitor, a ``tof.Detector``, or a chopper,
        a ``tof.Chopper``, the only components a tof run has a reading of.
    transmitted:
        True counts the neutrons a component lets through: of a chopper, those that pass its
        slits, as a McStas monitor just after it would see them. False counts all that reach
        it, also those a chopper blocks. A monitor blocks none: both count the same.

    Returns
    -------
    :
        The neutrons per bin of each component, by name, in counts, with their Poisson
        variances: N for a bin of N neutrons, an error of √N. Each has one dimension per entry of
        its bins, in their order, and its bins as its coordinates, in the units given.
    """
    return {name: _read_tof_map(name, tof_run, edges, transmitted) for name, edges in bins.items()}


def _read_tof_map(
    name: str, tof_run: tof.Result, bins: Mapping[str, sc.Variable], transmitted: bool
) -> sc.DataArray:
    """Count the neutrons of the first pulse at one component of a tof run."""
    neutrons = tof_run[name].data["pulse", 0]
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
