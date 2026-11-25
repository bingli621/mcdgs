"""The monitors of CSPEC.

The sample monitor is at ``sample_pos``, before the sample, and leaves the neutrons unchanged.
CSPEC_nl-specs_1.E.xlsx has no monitors: its size and bins are those of the sample monitor of
T-Rex.
"""

from mcdgs.core.components import Monitor_nD
from mcdgs.core.components.monitors import time_lambda_options

from .guides.straight import SAMPLE_PATH

SAMPLE_MONITOR_WIDTH = 0.1  # m, FIXME, not in CSPEC_nl-specs_1.E.xlsx: wider than the beam
SAMPLE_MONITOR_HEIGHT = 0.1  # m, FIXME, not in CSPEC_nl-specs_1.E.xlsx: wider than the beam
SAMPLE_MONITOR_LAMBDA_RESOLUTION = 0.02  # Å, FIXME, as the sample monitor of T-Rex
SAMPLE_MONITOR_TIME_RESOLUTION = 5e-5  # s, FIXME, as the sample monitor of T-Rex


def sample_monitor(
    time_limits: tuple[float, float] | None = None,
    lambda_limits: tuple[float, float] | None = None,
    lambda_resolution: float = SAMPLE_MONITOR_LAMBDA_RESOLUTION,
    time_resolution: float = SAMPLE_MONITOR_TIME_RESOLUTION,
) -> Monitor_nD:
    """Create the time and wavelength monitor at the sample position.

    It records the beam reaching the sample as one 2D map of time against wavelength, and leaves
    the neutrons unchanged. It is at ``sample_pos``, before the sample.

    Parameters
    ----------
    time_limits, lambda_limits, lambda_resolution, time_resolution:
        As in :func:`.time_lambda_options`, at ``SAMPLE_PATH``.

    Returns
    -------
    :
        A new monitor.
    """
    return Monitor_nD(
        xwidth=SAMPLE_MONITOR_WIDTH,
        yheight=SAMPLE_MONITOR_HEIGHT,
        restore_neutron=1,
        options=time_lambda_options(
            SAMPLE_PATH, time_limits, lambda_limits, lambda_resolution, time_resolution
        ),
        filename="sample_monitor",
        relative="sample_pos",
    )
