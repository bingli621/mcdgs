"""The detectors of T-Rex, after the sample.

In beam order, all on ``sample_pos``: the beamstop detector, a flat monitor on the beam axis
just in front of the banana detector; the beamstop, an absorbing plate of the same size at the
same place, which stops the direct beam there; and the banana detector, a cylinder around the
sample that lists every neutron reaching it.

The beamstop detector leaves the neutrons unchanged, so the beamstop behind it stops those it
counted. A neutron that misses both, or is scattered backward, goes on to the banana detector.

None of these is in Table of Optics v3.2. The banana detector takes its size and pixels from the
earlier McStas model of T-Rex.
"""

import math

from mcdgs.core.components import Beamstop, Monitor_nD

from .guides.straight import SAMPLE_Z
from .monitors import (
    PRE_SAMPLE_MONITOR_HEIGHT,
    PRE_SAMPLE_MONITOR_WIDTH,
    time_lambda_options,
)

BANANA_RADIUS = 3.0  # m, FIXME, not in Table of Optics v3.2: the sample-detector distance
BANANA_PHI_LIMITS = (-25, 15)  # deg, FIXME, not in Table of Optics v3.2: below, above the beam
BANANA_BOTTOM, BANANA_TOP = (
    BANANA_RADIUS * math.tan(math.radians(phi)) for phi in BANANA_PHI_LIMITS
)
BANANA_HEIGHT = BANANA_TOP - BANANA_BOTTOM  # m, 2.2028
BANANA_CENTER_Y = (BANANA_TOP + BANANA_BOTTOM) / 2  # m, -0.2975: below the beam
BANANA_THETA_LIMITS = (-36, 144)  # deg, FIXME, not in Table of Optics v3.2
BANANA_THETA_BINS = 288  # FIXME, not in Table of Optics v3.2: 0.625° pixels
BANANA_Y_BINS = 136  # FIXME, not in Table of Optics v3.2: pixels about 16 mm high

BEAMSTOP_GAP = 0.001  # m, FIXME, not in Table of Optics v3.2: in front of the banana detector
BEAMSTOP_DISTANCE = BANANA_RADIUS - BEAMSTOP_GAP  # m, from the sample
BEAMSTOP_WIDTH = PRE_SAMPLE_MONITOR_WIDTH  # m, FIXME, as the pre-sample monitor
BEAMSTOP_HEIGHT = PRE_SAMPLE_MONITOR_HEIGHT  # m, FIXME, as the pre-sample monitor
BEAMSTOP_LAMBDA_RESOLUTION = 0.02  # Å, FIXME, as the beam monitor
BEAMSTOP_TIME_RESOLUTION = 5e-5  # s, FIXME, as the beam monitor


def beamstop_detector(
    time_limits: tuple[float, float] | None = None,
    lambda_limits: tuple[float, float] | None = None,
    lambda_resolution: float = BEAMSTOP_LAMBDA_RESOLUTION,
    time_resolution: float = BEAMSTOP_TIME_RESOLUTION,
) -> Monitor_nD:
    """Create the time and wavelength monitor on the face of the beamstop.

    It records the beam the beamstop stops, mostly that through the sample, as one 2D map of
    time against wavelength, and leaves the neutrons unchanged. It is ``BEAMSTOP_DISTANCE``
    after ``sample_pos``, on the beam axis.

    Parameters
    ----------
    time_limits, lambda_limits, lambda_resolution, time_resolution:
        As in :func:`.time_lambda_options`, at ``SAMPLE_Z + BEAMSTOP_DISTANCE`` from the
        moderator.

    Returns
    -------
    :
        A new monitor.
    """
    return Monitor_nD(
        xwidth=BEAMSTOP_WIDTH,
        yheight=BEAMSTOP_HEIGHT,
        restore_neutron=1,
        options=time_lambda_options(
            SAMPLE_Z + BEAMSTOP_DISTANCE,
            time_limits,
            lambda_limits,
            lambda_resolution,
            time_resolution,
        ),
        filename="beamstop_detector",
        at=[0, 0, BEAMSTOP_DISTANCE],
        relative="sample_pos",
    )


def beamstop() -> Beamstop:
    """Create the beamstop, at the place and of the size of the beamstop detector.

    Returns
    -------
    :
        A new beamstop.
    """
    return Beamstop(
        xwidth=BEAMSTOP_WIDTH,
        yheight=BEAMSTOP_HEIGHT,
        at=[0, 0, BEAMSTOP_DISTANCE],
        relative="sample_pos",
    )


def banana_detector() -> Monitor_nD:
    """Create the banana detector, a cylinder of ``BANANA_RADIUS`` around the sample.

    Its axis is vertical, through ``sample_pos``; it spans ``BANANA_THETA_LIMITS`` in the
    horizontal plane, 0° along the beam. It is lowered by ``-BANANA_CENTER_Y`` so that, seen from
    the sample, it spans ``BANANA_PHI_LIMITS`` vertically, 0° on the beam. Its pixels
    are ``BANANA_THETA_BINS`` in angle by ``BANANA_Y_BINS`` in height, numbered from 0 as Mantid
    does. It lists every neutron that reaches it, with its pixel, time and wavelength, for
    :func:`mcdgs.mcstas.nexus_io.read_pixels`.

    Returns
    -------
    :
        A new monitor.
    """
    theta_start, theta_stop = BANANA_THETA_LIMITS
    return Monitor_nD(
        radius=BANANA_RADIUS,
        yheight=BANANA_HEIGHT,
        restore_neutron=1,
        options=(  # theta in deg: the pixels, then what each event lists
            f"mantid banana, theta limits=[{theta_start} {theta_stop}] bins={BANANA_THETA_BINS}, "
            f"y bins={BANANA_Y_BINS}, neutron pixel min=0 t lambda, list all neutrons"
        ),
        filename="banana_detector",
        at=[0, BANANA_CENTER_Y, 0],
        relative="sample_pos",
    )
