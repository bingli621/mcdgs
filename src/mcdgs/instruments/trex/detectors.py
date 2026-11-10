"""The detector of T-Rex, after the sample.

The banana detector, a cylinder around the sample on ``sample_pos`` that lists every neutron
reaching it. A neutron that misses the beamstop of :mod:`.monitors`, or is scattered backward,
reaches it.

It is not in Table of Optics v3.2: it takes its size and pixels from the earlier McStas model of
T-Rex.
"""

from mcdgs.core.components import Monitor_nD

BANANA_RADIUS = 3.0  # m, FIXME, not in Table of Optics v3.2: the sample-detector distance
BANANA_HEIGHT = 2.2  # m, FIXME, not in Table of Optics v3.2: -25 to 15 deg at 3 m
BANANA_CENTER_Y = -0.2978  # m, FIXME, not in Table of Optics v3.2: below the beam
BANANA_THETA_LIMITS = (-36, 144)  # deg, FIXME, not in Table of Optics v3.2
BANANA_THETA_BINS = 288  # FIXME, not in Table of Optics v3.2: 0.625° pixels
BANANA_Y_BINS = 136  # FIXME, not in Table of Optics v3.2: pixels about 16 mm high


def banana_detector() -> Monitor_nD:
    """Create the banana detector, a cylinder of ``BANANA_RADIUS`` around the sample.

    Its axis is vertical, through ``sample_pos``; it spans ``BANANA_THETA_LIMITS`` in the
    horizontal plane, 0° along the beam. It is ``BANANA_HEIGHT`` high and lowered by
    ``-BANANA_CENTER_Y``, so that, seen from the sample, it spans 25° below the beam to 15° above.
    Its pixels are ``BANANA_THETA_BINS`` in angle by ``BANANA_Y_BINS`` in height, numbered from 0
    as Mantid does. It lists every neutron that reaches it, with its pixel, time and wavelength, for
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
