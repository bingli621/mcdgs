"""The neutron source."""

from dataclasses import dataclass
from typing import Any

from .base import Component


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
