"""The detector of CSPEC, after the sample.

Twelve flat panels on a circle of ``PANEL_DISTANCE`` around the sample on ``sample_pos``, each
listing every neutron reaching it. They are the twelve ``Detector_tubes`` modules of drcspec, each
a pack of 32 vertical tubes of 350 pixels, made ``Monitor_nD`` rectangles of the same size, place
and pixels.

CSPEC_nl-specs_1.E.xlsx has no detector: the panels take their size, places and pixels from
drcspec.
"""

import math

from mcdgs.core.components import Monitor_nD

PANEL_DISTANCE = 3.5  # m, FIXME, not in CSPEC_nl-specs_1.E.xlsx: det_rad, sample to panel center
PANEL_WIDTH = 0.823239966741101  # m, FIXME, not in CSPEC_nl-specs_1.E.xlsx: 32 tubes side by side
PANEL_HEIGHT = 3.5  # m, FIXME, not in CSPEC_nl-specs_1.E.xlsx: det_height, the tube length
PANEL_X_BINS = 32  # FIXME, not in CSPEC_nl-specs_1.E.xlsx: the tubes of a module, N
PANEL_Y_BINS = 350  # FIXME, not in CSPEC_nl-specs_1.E.xlsx: the pixels of a tube, no, 1 cm
# FIXME, not in CSPEC_nl-specs_1.E.xlsx: the rotation about +y of the arm of each module, in deg;
# the module is on the +x axis of its arm, so its center is 90 deg further, from the beam
PANEL_ARM_ANGLES = (
    47.11249999999998,
    33.067499999999995,
    19.022499999999994,
    4.977499999999992,
    -9.067499999999995,
    -23.11250000000001,
    -37.1575,
    -51.2025,
    -65.2475,
    -79.2925,
    -100.7075,
    -114.75249999999998,
)
PANEL_PIXELS = PANEL_X_BINS * PANEL_Y_BINS  # pixels of a panel; panel i starts at i times this


def panel_angle(index: int) -> float:
    """Return the angle of the center of a panel from the beam, about +y.

    Parameters
    ----------
    index:
        The panel, from 0 to 11, from the largest angle to the smallest.

    Returns
    -------
    :
        The angle in deg, 0 along the beam, positive toward +x, that of its arm plus 90.
    """
    return PANEL_ARM_ANGLES[index] + 90


def detector_panel(index: int) -> Monitor_nD:
    """Create a detector panel, a flat rectangle facing the sample.

    Its center is ``PANEL_DISTANCE`` from ``sample_pos`` in the horizontal plane, at
    :func:`panel_angle` from the beam, and it faces the sample. It is ``PANEL_WIDTH`` wide and
    ``PANEL_HEIGHT`` high, centered on the height of the beam. Its pixels are ``PANEL_X_BINS``
    tubes across by ``PANEL_Y_BINS`` along each tube, numbered as Mantid does from
    ``index * PANEL_PIXELS``, so that the pixels of all panels differ. It lists every neutron that
    reaches it, with its pixel, time and wavelength, for :func:`mcdgs.mcstas.reader.read_pixels`.

    Parameters
    ----------
    index:
        The panel, from 0 to 11, from the largest angle to the smallest.

    Returns
    -------
    :
        A new monitor.
    """
    angle = panel_angle(index)
    phi = math.radians(angle)
    return Monitor_nD(
        xwidth=PANEL_WIDTH,
        yheight=PANEL_HEIGHT,
        restore_neutron=1,
        options=(  # the pixels, then what each event lists
            f"mantid square, x bins={PANEL_X_BINS}, y bins={PANEL_Y_BINS}, "
            f"neutron pixel min={index * PANEL_PIXELS} t lambda, list all neutrons"
        ),
        filename=f"detector_panel_{index}",
        at=[PANEL_DISTANCE * math.sin(phi), 0, PANEL_DISTANCE * math.cos(phi)],
        rotated=[0, angle, 0],
        relative="sample_pos",
    )


def detector_panels() -> dict[str, Monitor_nD]:
    """Create the twelve detector panels of :func:`detector_panel`.

    Returns
    -------
    :
        Each panel by name, ``detector_panel_0`` to ``detector_panel_11``, from the largest angle
        to the smallest.
    """
    return {f"detector_panel_{i}": detector_panel(i) for i in range(len(PANEL_ARM_ANGLES))}
