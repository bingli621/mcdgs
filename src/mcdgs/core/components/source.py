"""The neutron sources: ESS_butterfly and Polygon_ESS_butterfly.

``Polygon_ESS_butterfly`` is ``ESS_butterfly`` emitting only where a chopper train lets neutrons
through.
"""

from dataclasses import dataclass, field
from typing import Any

from .base import Component
from .choppers import ChopperTrainParameters

ESS_SOURCE_FREQUENCY = 14  # Hz, the pulse frequency
ESS_SOURCE_LENGTH = 3.5e-3  # s, the pulse length
# the widest band the source simulates: the band is clipped to it
L_MIN = 0.1  # Å
L_MAX = 10.0  # Å


@dataclass(kw_only=True, repr=False)  # repr: Component's
class ESS_butterfly(Component):
    """The ESS butterfly moderator, as McStas ``ESS_butterfly``.

    Name and instance keywords (``name``, ``at``, ...) are as in :class:`.Component`.

    Parameters
    ----------
    Lmin, Lmax:
        Wavelength band simulated in Å, or the instrument parameters that provide it.
    sector:
        Sector of the beamline, "N", "S", "E" or "W".
    beamline:
        Beamline number in the sector, 1 to 10 or 11 depending on the sector.
    acc_power:
        Accelerator power in MW.
    yheight:
        Moderator height in m, 0.03 or 0.06.
    cold_frac:
        Fraction of the events emitted from the cold moderator, 0 to 1.
    dist:
        Distance in m from the source to the focusing rectangle, at (0, 0, ``dist``).
    focus_xw, focus_yh:
        Width and height of the focusing rectangle in m.
    n_pulses:
        Number of pulses simulated; 0 and 1 are one pulse.
    c_performance, t_performance:
        Factors on the cold and the thermal brilliance, > 0.
    target_index:
        Index of the component to focus on, relative to this one (+1 is the next); McStas then
        computes ``dist`` from it. 0 uses ``dist``.
    tmax_multiplier:
        Maximum emission time, in units of the ESS pulse duration.
    tfocus_dist:
        Position in m of the time-focusing window along z. 0 is no time focusing.
    tfocus_time, tfocus_width:
        Time and width of the time-focusing window in s.

    Notes
    -----
    ``Lmin`` and ``Lmax`` are required, as in McStas. ``sector``, ``beamline`` and
    ``acc_power`` default to ``"W"``, 7 and 2 MW; every other default is McStas's.

    Examples
    --------
    >>> source = ESS_butterfly(
    ...     name="source",
    ...     sector="W",
    ...     beamline=7,
    ...     acc_power=2,
    ...     yheight=0.03,
    ...     cold_frac=0.3,
    ...     dist=2,
    ...     focus_xw=0.095,
    ...     focus_yh=0.035,
    ...     n_pulses=1,
    ...     Lmin=0.1,
    ...     Lmax=10.0,
    ... )
    >>> source.mcstas_parameters()["sector"]
    'W'
    """

    Lmin: float | str
    Lmax: float | str
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

    @property
    def frequency(self) -> float:
        """Pulse frequency in Hz, for an instrument's choppers to be phased to.

        McStas fixes it at 14 Hz (``ESS_SOURCE_FREQUENCY`` in ess_source-lib.h), so it is read
        only, and no McStas parameter.
        """
        return ESS_SOURCE_FREQUENCY

    @property
    def pulse_length(self) -> float:
        """Pulse length in s, for an instrument's choppers to be phased to.

        3.5 ms, ``ESS_SOURCE_LENGTH``. It is read only, and no McStas parameter: McStas's own
        pulse, ``ESS_SOURCE_DURATION`` in ess_source-lib.h, is 2.857 ms.
        """
        return ESS_SOURCE_LENGTH

    def mcstas_parameters(self) -> dict[str, Any]:
        """Return the McStas parameters, once ``sector`` is checked.

        Returns
        -------
        :
            The fields.

        Raises
        ------
        ValueError
            If ``sector`` is not N, S, E or W.
        """
        if self.sector not in ("N", "S", "E", "W"):
            raise ValueError(f"{self.name}: sector {self.sector!r} is not 'N', 'S', 'E' or 'W'")
        return super().mcstas_parameters()


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
        As in :class:`ESS_butterfly`.
    choppers:
        The disks of the chopper train, in any order, as :class:`.ChopperTrainParameters`.
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
    and ``acc_power`` default to ``"W"``, 7 and 2 MW, as in :class:`ESS_butterfly`; every other
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
        """Pulse frequency in Hz, as :attr:`ESS_butterfly.frequency`."""
        return ESS_SOURCE_FREQUENCY

    @property
    def pulse_length(self) -> float:
        """Pulse length in s, as :attr:`ESS_butterfly.pulse_length`."""
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
