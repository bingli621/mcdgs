"""The incident energy of neutrons from a straight flight between two points.

``EI_GRAPH`` is a scipp graph built from the kernels of scippneutron. It takes a straight flight
from a start to an end, such as from one monitor to another, or from the source to a monitor,
and the time at each, and gives the flight path, the time of flight, the incident wavelength and
the incident energy. :func:`incident_flight` follows it for given positions and times.
"""

from typing import cast

import scipp as sc
from scipp.coords.graph import GraphDict
from scipp.coords.rule import Kernel
from scippneutron.conversion import beamline as beamline_kernels
from scippneutron.conversion import tof as tof_kernels

__all__ = ["EI_GRAPH", "incident_flight", "time_of_flight"]


# --8<-- [start:ei_graph]
def time_of_flight(time_at_start: sc.Variable, time_at_end: sc.Variable) -> sc.Variable:
    """Return the time of flight from the start to the end of a straight flight."""
    return time_at_end - time_at_start


# scippneutron types it for data arrays and datasets too, scipp's graph only for variables
_flight_path = cast(Kernel, beamline_kernels.total_straight_beam_length_no_scatter)

# a straight flight from a start to an end, as scippneutron's beam from a source to a position,
# without scattering
EI_GRAPH: GraphDict = {
    "Ltotal": _flight_path,  # from source_position to position
    "source_position": "start_position",
    "position": "end_position",
    "tof": time_of_flight,
    "incident_wavelength": tof_kernels.wavelength_from_tof,  # from tof and Ltotal
    "incident_energy": tof_kernels.incident_energy_from_wavelength,
}
# --8<-- [end:ei_graph]


def incident_flight(
    start_position: sc.Variable,
    time_at_start: sc.Variable,
    end_position: sc.Variable,
    time_at_end: sc.Variable,
) -> dict[str, sc.Variable]:
    """Compute the incident energy of a straight flight with ``EI_GRAPH``.

    Parameters
    ----------
    start_position, end_position:
        Where the flight starts and ends, as vectors.
    time_at_start, time_at_end:
        When it passes each. The times may have a dimension, such as one time per sub-pulse,
        and the time at the start may be one for all.

    Returns
    -------
    :
        ``Ltotal``, ``tof``, ``incident_wavelength`` and ``incident_energy``, by name, with the
        dimensions of the times.

    Examples
    --------
    A neutron of 2 Å flies at 1978 m/s, so 10 m take 5.0556 ms:

    >>> flight = incident_flight(
    ...     sc.vector([0.0, 0.0, 0.0], unit="m"),
    ...     sc.scalar(0.0, unit="ms"),
    ...     sc.vector([0.0, 0.0, 10.0], unit="m"),
    ...     sc.scalar(5.05557, unit="ms"),
    ... )
    >>> print(f"{flight['incident_wavelength'].to(unit='Å').value:.4f} Å")
    2.0000 Å
    >>> print(f"{flight['incident_energy'].to(unit='meV').value:.2f} meV")
    20.45 meV
    """
    sizes = (time_at_start + time_at_end).sizes  # either may have the dimension
    flight = sc.DataArray(
        sc.zeros(sizes=sizes),  # only the coordinates count
        coords={
            "start_position": start_position,
            "time_at_start": time_at_start,
            "end_position": end_position,
            "time_at_end": time_at_end,
        },
    )
    names = ["Ltotal", "tof", "incident_wavelength", "incident_energy"]
    flight = cast(
        sc.DataArray, flight.transform_coords(names, graph=EI_GRAPH, keep_intermediate=False)
    )
    return {name: flight.coords[name] for name in names}
