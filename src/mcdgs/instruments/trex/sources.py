"""The source of T-Rex."""

from mcdgs.core.components import (
    ESS_butterfly,
)


def source() -> ESS_butterfly:
    """Create the source component for T-Rex.

    Returns
    -------
    :
        A new source.
    """
    return ESS_butterfly(
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
