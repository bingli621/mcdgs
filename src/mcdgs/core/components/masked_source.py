"""The masked source: Polygon_ESS_butterfly, and the disks of its chopper train.

``Polygon_ESS_butterfly`` is ``ESS_butterfly`` emitting only where a chopper train lets neutrons
through. :class:`ChopperTrainParameters` describes one disk of the train, as chopper-lib takes it.
"""

from dataclasses import dataclass, field
from typing import Any, TypedDict

from .base import Component
from .source import ESS_SOURCE_FREQUENCY, ESS_SOURCE_LENGTH


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
class Polygon_ESS_butterfly(Component):
    """The ESS butterfly source masked by a chopper train, as McStas ``Polygon_ESS_butterfly``.

    Before the run, it computes the exact region of (inverse velocity, emission time) that the
    disks of ``choppers`` let through, a union of convex polygons. It then emits only inside
    it: with ``resample``, a ray drawn outside is drawn again inside, and every ray's weight is
    multiplied by the share of the sampled plane the region covers.

    Name and instance keywords (``name``, ``at``, ...) are as in :class:`.Component`.

    Parameters
    ----------
    Lmin, Lmax, sector, beamline, acc_power, yheight, cold_frac, dist, focus_xw, focus_yh, \
n_pulses, c_performance, t_performance, target_index, tmax_multiplier, tfocus_dist, \
tfocus_time, tfocus_width:
        As in :class:`.ESS_butterfly`.
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
    ``Lmin`` and ``Lmax`` are required, as in McStas. ``sector``, ``beamline``
    and ``acc_power`` default to ``"W"``, 7 and 2 MW, as in :class:`.ESS_butterfly`; every other
    default is McStas's.

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

    Lmin: float | str
    Lmax: float | str
    choppers: list[ChopperTrainParameters] = field(default_factory=list)
    sector: str = "W"
    beamline: int = 7
    acc_power: float = 2
    yheight: float = 0.03
    cold_frac: float = 0.5
    dist: float = 0
    focus_xw: float = 0
    focus_yh: float = 0
    n_pulses: int = 1
    c_performance: float = 1
    t_performance: float = 1
    target_index: int = 0
    tmax_multiplier: float = 3
    tfocus_dist: float = 0
    tfocus_time: float = 0
    tfocus_width: float = 0
    filename: str | None = None
    path_spread_fraction: float = 0
    noise_fraction: float = 0
    emission_depth: int = 1
    use_region: int = 1
    resample: int = 1
    save_polygons: int = 1
    verify_acceptance: int = 1

    @property
    def frequency(self) -> float:
        """Pulse frequency in Hz, as :attr:`.ESS_butterfly.frequency`."""
        return ESS_SOURCE_FREQUENCY

    @property
    def pulse_length(self) -> float:
        """Pulse length in s, as :attr:`.ESS_butterfly.pulse_length`."""
        return ESS_SOURCE_LENGTH

    def mcstas_parameters(self) -> dict[str, Any]:
        """Return the McStas parameters, with ``resample`` 0 where McStas refuses 1.

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
