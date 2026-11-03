"""The source of T-Rex."""

from dataclasses import dataclass

from mcdgs.core.components import (
    ESS_butterfly,
)


@dataclass(kw_only=True, repr=False)  # repr: Component's
class TRexSource(ESS_butterfly):
    """The ESS butterfly moderator of T-Rex, with the pulse frequency its choppers follow.

    Parameters are those of :class:`.ESS_butterfly`, plus:

    Parameters
    ----------
    frequency:
        Pulse frequency in Hz, for the choppers to be phased to. Not a McStas parameter: McStas
        fixes it at 14 Hz (``ESS_SOURCE_FREQUENCY``), so the builder leaves it out and changing
        it does not change the pulses.

    Notes
    -----
    The McStas component is still ``ESS_butterfly``, the class this one extends.
    """

    frequency: float = 14  # not for McStas, which fixes it at this value


def source() -> TRexSource:
    """Create the source component for T-Rex.

    Returns
    -------
    :
        A new source.
    """
    return TRexSource(
        sector="W",
        beamline=7,
        yheight=0.03,
        cold_frac=0.3,
        dist=2,
        focus_xw=0.095,
        focus_yh=0.035,
        Lmin=0.1,
        Lmax=10.0,
        n_pulses=1,
        acc_power=2,
    )
