"""Read the data of the components of a McStas run, as McStasScript loads them.

:func:`read_mcstas_map` reads the time and wavelength map of a component, in n/s per bin, on its
bin edges (:func:`.monitor_bins`), as a scipp data array with its errors.
:func:`.read_tof_map` reads the same map of a tof run.

The neutrons of a detector are read from the NeXus file instead: see :mod:`.nexus_io`.
"""

from collections.abc import Mapping

import mcstasscript as ms
import numpy as np
import scipp as sc
from mcstasscript.data.data import McStasData


def read_mcstas_map(
    name: str, mcstas_data: list[McStasData], bins: Mapping[str, sc.Variable]
) -> sc.DataArray:
    """Read the time and wavelength map of a component of one McStas run, in n/s per bin.

    Parameters
    ----------
    name:
        The name of the component, such as a monitor.
    mcstas_data:
        The data of the components of the McStas run, as ``mcstasscript.load_data`` reads them.
        That of this one must be a 2D map of time against wavelength,
        ``"time ... lambda ..."``, on ``bins``.
    bins:
        The bin edges of the map, ``"toa"`` and ``"wavelength"``, in that order.

    Returns
    -------
    :
        The intensity per bin, in n/s, with the square of its error as variances. Its
        dimensions are those of ``bins``, in their order, (time bins, wavelength bins), and its
        coordinates the bin edges.

    Raises
    ------
    ValueError
        If several data sets have the name, or the map does not have the shape of ``bins``.
    """
    component = ms.name_search(name, mcstas_data)
    if isinstance(component, list):  # McStasScript returns all of them, as a list
        # the name is ambiguous, not of a wrong type: a ValueError
        raise ValueError(f"{name}: {len(component)} data sets of that name, not one")  # noqa: TRY004
    # McStasScript has the map as (wavelength, time), the y axis first
    intensity = np.asarray(component.Intensity, dtype=float).T
    error = np.asarray(component.Error, dtype=float).T
    shape = tuple(len(edges) - 1 for edges in bins.values())
    if intensity.shape != shape:
        raise ValueError(
            f"{name}: McStas has a map of {intensity.shape} bins, the bins are {shape}"
        )
    return sc.DataArray(
        sc.array(dims=list(bins), values=intensity, variances=error**2, unit="1/s"),
        coords=dict(bins),
    )
