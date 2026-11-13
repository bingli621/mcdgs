"""The masked source: Polygon_ESS_butterfly, and the disks of its chopper train.

``Polygon_ESS_butterfly`` is ``ESS_butterfly`` emitting only where a chopper train lets neutrons
through. :class:`ChopperTrainParameters` describes one disk of the train, as chopper-lib takes it,
and :func:`chopper_train` describes the disks an instrument passes it, each with its flight path.
"""

import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any, TypedDict

from .choppers import NXdisk_chopper
from .source import ESS_butterfly


class ChopperTrainParameters(TypedDict):
    """One disk of the chopper train of a :class:`Polygon_ESS_butterfly`.

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


@dataclass(kw_only=True, repr=False)  # repr: Component's
class Polygon_ESS_butterfly(ESS_butterfly):
    """The ESS butterfly source masked by a chopper train, as McStas ``Polygon_ESS_butterfly``.

    It extends :class:`.ESS_butterfly` as the McStas component extends it: it has all of its
    parameters, defaults, ``frequency`` and ``pulse_length``, and checks ``sector`` the same way.

    Before the run, it computes the exact region of (inverse velocity, emission time) that the
    disks of ``choppers`` let through, a union of convex polygons. It then emits only inside
    it: with ``resample``, a ray drawn outside is drawn again inside, and every ray's weight is
    multiplied by the share of the sampled plane the region covers.

    Name and instance keywords (``name``, ``at``, ...) are as in :class:`.Component`.

    Parameters
    ----------
    choppers:
        The disks of the chopper train, in any order, as :class:`ChopperTrainParameters`.
        McStas's ``chopper_count`` is computed from it. Without a disk, the region is all the
        source samples, so the source is ``ESS_butterfly``.
    filename:
        Base name of the output file. If ``None``, the name of the component.
    path_spread_fraction:
        Extra flight path a ray may have, as a fraction of each disk's ``path``, >= 0.
    noise_fraction:
        Probability of keeping a ray outside the region, 0 to 1.
    emission_depth:
        1 allows for the depth of the moderator surfaces along the beam, 0 emits as though from
        the component's origin.
    use_region:
        1 limits the emission to the region, 0 does not.
    resample:
        1 draws a ray outside the region again, inside it; 0 absorbs it. With time focusing,
        ``use_region`` 0 or a ``noise_fraction``, McStas refuses 1, so the source runs with 0
        and says so.
    save_polygons:
        1 writes the region to ``{filename}.json``.
    verify_acceptance:
        1 counts the rays drawn and those inside the region, and reports their ratio next to the
        acceptance the geometry gives.

    Notes
    -----
    Every default of its own is McStas's.

    Examples
    --------
    >>> source = Polygon_ESS_butterfly(
    ...     name="source",
    ...     Lmin=1.0,
    ...     Lmax=3.0,
    ...     choppers=[
    ...         {
    ...             "speed": "nu_ps1",
    ...             "delay": "delay_ps1",
    ...             "beam": 180,
    ...             "edges": [-10, 10],
    ...             "path": 107.95,
    ...             "aperture": 3.2,
    ...         }
    ...     ],
    ... )
    >>> source.mcstas_parameters()["chopper_count"]
    1
    """

    choppers: list[ChopperTrainParameters] = field(default_factory=list)
    filename: str | None = None
    path_spread_fraction: float = 0
    noise_fraction: float = 0
    emission_depth: int = 1
    use_region: int = 1
    resample: int = 1
    save_polygons: int = 1
    verify_acceptance: int = 1

    def mcstas_parameters(self) -> dict[str, Any]:
        """Return the McStas parameters, with ``resample`` 0 where McStas refuses 1.

        ``sector`` is checked as in :meth:`.ESS_butterfly.mcstas_parameters`.

        McStas stops the run in INITIALIZE, after the compilation, if ``resample`` is 1 with
        time focusing, ``use_region`` 0 or a ``noise_fraction``. Then ``resample`` is 0 here,
        and the field stays as it is, and a line printed says why.

        An unset ``filename`` is the component's name: McStas names the file after the
        component only for an empty one, and McStas's ``"NULL"`` would name it ``NULL.json``.

        Returns
        -------
        :
            The fields, with ``filename`` set and ``resample`` as McStas takes it, plus
            ``chopper_count``.

        Raises
        ------
        ValueError
            If ``sector`` is not N, S, E or W.
        """
        parameters = super().mcstas_parameters() | {
            "chopper_count": len(self.choppers),
            "filename": self.name if self.filename is None else self.filename,
        }
        refused = [
            reason
            for reason, applies in (
                ("time focusing", self.tfocus_width > 0),
                ("use_region=0", not self.use_region),
                (f"noise_fraction={self.noise_fraction}", self.noise_fraction != 0),
            )
            if applies
        ]
        if self.resample and refused:
            print(
                f"{self.name}: resample=1 does not go with {', '.join(refused)}; running with "
                "resample=0, which absorbs the rays outside the region"
            )
            parameters["resample"] = 0
        return parameters


def aperture(disk: NXdisk_chopper) -> float:
    """Return the angular width of the beam on a disk, about its spindle, as chopper-lib takes it.

    The inner corners of the beam window, nearest the spindle, are the furthest round the disk.
    NXdisk_chopper hangs the window from the height of the rim at its edges,
    ``sqrt(radius**2 - (xwidth/2)**2)``, so those corners are ``yheight`` inside it.

    Parameters
    ----------
    disk:
        The disk, with its beam window.

    Returns
    -------
    :
        ``2 atan2(xwidth / 2, sqrt(radius**2 - (xwidth / 2)**2) - yheight)``, in deg.

    Examples
    --------
    >>> disk = NXdisk_chopper(slit_edges=[-10, 10], radius=0.5, yheight=0.1, xwidth=0.1)
    >>> round(aperture(disk), 3)
    14.339
    """
    half_width = disk.xwidth / 2
    reach = math.sqrt(disk.radius**2 - half_width**2) - disk.yheight
    return math.degrees(2 * math.atan2(half_width, reach))


def chopper_train(disks: Sequence[tuple[NXdisk_chopper, float]]) -> list[ChopperTrainParameters]:
    """Describe disks to chopper-lib, as :class:`Polygon_ESS_butterfly` takes them.

    Each entry has the disk's edges, beam angle and frequency and delay, which may be instrument
    parameters, and its :func:`aperture`. Both McStas components put the disk angle ``a`` on the
    beam axis at ``delay + (beam_angle - a) / (360 nu)``, so the delay is the same.

    Parameters
    ----------
    disks:
        Each disk with its flight path from the moderator, in m: the disk is placed relative to
        other components, so the instrument gives the path.

    Returns
    -------
    :
        One entry per disk, in the order of ``disks``.

    Raises
    ------
    TypeError
        If a disk's ``beam_angle`` is an instrument parameter: chopper-lib takes a number.

    Examples
    --------
    >>> disk = NXdisk_chopper(slit_edges=[-10, 10], nu="nu_ps1", delay="delay_ps1", beam_angle=180)
    >>> chopper_train([(disk, 107.95)])[0]["beam"], chopper_train([(disk, 107.95)])[0]["path"]
    (180, 107.95)
    """
    train: list[ChopperTrainParameters] = []
    for disk, path in disks:
        if isinstance(disk.beam_angle, str):
            raise TypeError(f"beam_angle {disk.beam_angle!r} is not a number")
        train.append(
            {
                "speed": disk.nu,
                "delay": disk.delay,
                "beam": disk.beam_angle,
                "edges": list(disk.slit_edges),
                "path": path,
                "aperture": aperture(disk),
            }
        )
    return train
