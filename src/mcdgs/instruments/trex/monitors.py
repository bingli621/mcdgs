"""The monitors of T-Rex.

Each monitor leaves the neutrons unchanged. The beam monitor sits in the middle of the bunker gap
of row 29, after guide 4, on the arm of guide 4; the sample monitor at ``sample_pos``.
"""

from mcdgs.core.components import Monitor_nD

from .guides.bispectral import GUIDE_4_TILT, GUIDE_4A_IN_Z
from .guides.extraction import z_along_arm

BEAM_MONITOR_GAP_Z = 5.892  # m, cell B29
BEAM_MONITOR_GAP_LENGTH = 0.02  # m, cell D29
BEAM_MONITOR_Z = BEAM_MONITOR_GAP_Z + BEAM_MONITOR_GAP_LENGTH / 2  # m, mid-gap


def beam_monitor(
    time_limits: tuple[float, float, int], lambda_limits: tuple[float, float, int]
) -> Monitor_nD:
    """Create the time and wavelength monitor in the middle of the bunker gap of row 29.

    It records the intensity leaving guide 4 as one 2D map of time against wavelength, and
    leaves the neutrons unchanged.

    Parameters
    ----------
    time_limits:
        The time axis: start and stop in s, and the number of bins.
    lambda_limits:
        The wavelength axis: start and stop in Å, and the number of bins.

    Returns
    -------
    :
        A new monitor.
    """
    t_start, t_stop, t_bins = time_limits
    l_start, l_stop, l_bins = lambda_limits

    return Monitor_nD(
        xwidth=0.06,  # FIXME cell E32, the entrance of guide 5
        yheight=0.046260695185495173,  # FIXME cell F32, the entrance of guide 5
        restore_neutron=1,
        options=(  # time in s, lambda in Å: one 2D map
            f"time limits=[{t_start} {t_stop}] bins={t_bins} "
            f"lambda limits=[{l_start} {l_stop}] bins={l_bins}"
        ),
        filename="beam_monitor",
        at=[0, 0, z_along_arm(BEAM_MONITOR_Z, GUIDE_4A_IN_Z, GUIDE_4_TILT)],
        relative="arm_guide_4",
    )


def sample_monitor() -> Monitor_nD:
    """Create the wavelength monitor at the sample position.

    It records the beam reaching the sample and leaves the neutrons unchanged. It is at
    ``sample_pos``, before the sample.

    Returns
    -------
    :
        A new monitor.
    """
    return Monitor_nD(
        xwidth=0.1,  # FIXME, not in Table of Optics v3.2: wider than the beam
        yheight=0.1,  # FIXME, not in Table of Optics v3.2: wider than the beam
        restore_neutron=1,
        options="lambda limits=[0.5 5.3] bins=12",  # FIXME, not in Table of Optics v3.2
        filename="sample_monitor",
        relative="sample_pos",
    )
