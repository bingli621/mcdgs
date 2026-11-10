"""Choppers: NXdisk_chopper, and the disks of a chopper train as chopper-lib takes them."""

from dataclasses import dataclass
from typing import Any, TypedDict

from .base import Component


@dataclass(kw_only=True, repr=False)  # repr: Component's
class NXdisk_chopper(Component):
    """A disk chopper with arbitrary openings, as McStas ``NXdisk_chopper``.

    Name and instance keywords (``name``, ``at``, ...) are as in :class:`.Component`.

    Parameters
    ----------
    slit_edges:
        Opening and closing edge of each slit from the disk's zero mark, in pairs, increasing,
        in deg. McStas's ``n_edges`` is computed from it.
    radius:
        Outer radius of the disk in m.
    yheight:
        Radial clearance of the openings across the whole window in m. 0 reaches the spindle.
    xwidth:
        Chord width of the beam window in m. 0 leaves the opening unrestricted.
    nu:
        Signed rotation frequency in Hz. 0 parks the disk at ``park_angle``.
    delay:
        Time at which the disk's zero mark is on the beam, in s.
    park_angle:
        Where the zero mark stands while parked, in deg. Only read when ``nu`` is 0.
    beam_angle:
        Where the beam crosses the disk from the zero mark, in the sense of ``slit_edges``, in
        deg. 0 puts the beam at the top of the disk, 180 hangs the disk above the beam.
    zero_angle:
        Angle from the component's +y axis to the top-dead-center pickup, in deg.
    jitter:
        Timing jitter applied to each ray, in s.
    abs_out:
        If positive, rays that miss the disk entirely are absorbed.
    verbose:
        Output level: -1 silent, 0 errors, 1 warnings, 2 info, 3 details.

    Notes
    -----
    ``slit_edges`` is required, as in McStas. Every default is McStas's.

    Examples
    --------
    >>> chopper = NXdisk_chopper(name="chopper", slit_edges=[-10, 10], nu="f", delay=0.1)
    >>> parameters = chopper.mcstas_parameters()
    >>> parameters["slit_edges"], parameters["n_edges"], parameters["nu"]
    ([-10, 10], 2, 'f')
    """

    slit_edges: list[float]
    radius: float = 0.35
    yheight: float = 0
    xwidth: float = 0
    nu: float | str = 0
    delay: float | str = 0
    park_angle: float | str = 0
    beam_angle: float | str = 0
    zero_angle: float = 0
    jitter: float = 0
    abs_out: int = 0
    verbose: int = 0

    def mcstas_parameters(self) -> dict[str, Any]:
        """Return the McStas parameters, with the number of edges.

        Returns
        -------
        :
            The fields, plus ``n_edges``.
        """
        return super().mcstas_parameters() | {"n_edges": len(self.slit_edges)}


class ChopperTrainParameters(TypedDict):
    """One disk of the chopper train of a :class:`.Polygon_ESS_butterfly`.

    It is chopper-lib's ``chopper_parameters``, as a plain dict, so that a parameters file holds
    it as it is.

    Attributes
    ----------
    speed:
        Signed rotation frequency in Hz, or the instrument parameter that provides it.
    delay:
        Time at which the disk point at ``beam`` is on the beam, in s, or the instrument
        parameter that provides it.
    beam:
        Where the beam crosses the disk from its zero mark, in deg.
    edges:
        Opening and closing edge of each slit from the disk's zero mark, in pairs, increasing,
        in deg.
    path:
        Flight path from the moderator to the disk, in m.
    aperture:
        Angular width of the beam on the disk, about its spindle, in deg. 0 is a pencil beam.
    """

    speed: float | str
    delay: float | str
    beam: float
    edges: list[float]
    path: float
    aperture: float
