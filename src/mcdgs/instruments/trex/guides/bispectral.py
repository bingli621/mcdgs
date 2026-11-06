"""Guide 4 of T-Rex, the bispectral switch: thermal guide 4a, or the bender and cold guide 4b.

Dimensions are from the Table of Optics V3.2 (the .xlsx file in ``reference/``), whose rows are
named in the comments. The values are the sheet's, not rounded.

Guide 4 tapers at 0.4267 deg (cell N25), not at the 0.4300 deg of guides 1 to 3, and the sheet
carries the exit width of guide 3 unchanged across the gap, so guide 4 has its own arm, which keeps
the -x wall on the line of :mod:`.extraction`. Guides 4a and 4b lie on one axis and share that arm.

Guides 4a (thermal) and 4b (cold, behind the bender) occupy the same slot and share the same exit
plane: they are alternative optics, chosen by their ``when`` condition, not parallel beamlines.
The bender's entrance is oriented toward the cold moderator and its exit tangent is parallel to
the guide axis, so nothing downstream of it is rotated.

After guide 4, before guide 5, the bunker gap of row 29 is left for the ESS beam monitor, which
is in :mod:`mcdgs.instruments.trex.monitors`.
"""

import math

from mcdgs.core.components import Arm, Guide_gravity, Pol_bender

from .extraction import entrance_x, half_taper_angle, z_along_arm

BENDER_LENGTH = 0.05  # m, cell B25 minus cell B24
BENDER_RADIUS = 7.2  # FIXME, not in Table of Optics v3.2, in m
BENDER_ROTATION = -math.degrees(math.atan(BENDER_LENGTH / BENDER_RADIUS))  # −0.39788 deg

GUIDE_4A_IN_Z = 5.40065  # m, cell B24, also the entrance of the bender
GUIDE_4A_IN_WIDTH = 0.063838  # m, cell E24
GUIDE_4A_LENGTH = 0.481  # m, cell D24
GUIDE_4_OUT_WIDTH = 0.060255817017560886  # m, cells G24 and G25
GUIDE_4_TILT = half_taper_angle(
    GUIDE_4A_IN_WIDTH, GUIDE_4_OUT_WIDTH, GUIDE_4A_LENGTH
)  # -0.2134 deg, half of cells N25 and N26


def arm_guide_4() -> Arm:
    """Create the arm that guides 4a and 4b and the bender share.

    Returns
    -------
    :
        A new arm, on the bisector of the walls of guide 4.
    """
    return Arm(
        at=[entrance_x(GUIDE_4A_IN_WIDTH, GUIDE_4_TILT), 0, GUIDE_4A_IN_Z],
        rotated=[0, GUIDE_4_TILT, 0],
    )


def guide_4a(when: str | None = None) -> Guide_gravity:
    """Create guide 4a, the thermal guide, in place of the bender and guide 4b.

    Parameters
    ----------
    when:
        McStas condition for the guide to act, such as ``"bender_in==0"``. None, the default,
        always.

    Returns
    -------
    :
        A new guide.
    """
    return Guide_gravity(  # row 24
        w1=GUIDE_4A_IN_WIDTH,
        h1=0.044593536722227093,  # cell F24
        w2=GUIDE_4_OUT_WIDTH,
        h2=0.046152195301135549,  # cell H24
        l=GUIDE_4A_LENGTH,
        mleft=1.5,  # cell O24
        mright=1.5,  # cell P24
        mtop=3.5,  # cell Q24
        mbottom=3.5,  # cell R24
        relative="arm_guide_4",
        when=when,
    )


def bender(when: str | None = None) -> Pol_bender:
    """Create the bender, 63.84 mm x 46 mm x 50 mm (WxHxL), with 0.15 mm Si channels.

    Parameters
    ----------
    when:
        McStas condition for the bender to act, such as ``"bender_in==1"``. None, the default,
        always.

    Returns
    -------
    :
        A new bender.
    """
    # Si channels, NiTi coating. The mirror values are R0, Qc, alpha, m and W
    mirror = (0.99, 0.0219, 3.02, 4.0, 0.003)  # FIXME, not in Table of Optics v3.2

    return Pol_bender(
        xwidth=0.06384,  # FIXME, not in Table of Optics v3.2
        yheight=0.046,  # FIXME, not in Table of Optics v3.2
        length=BENDER_LENGTH,
        radius=BENDER_RADIUS,
        G=0,  # as every Guide_gravity
        nslit=426,  # FIXME, not in Table of Optics v3.2
        d=1e-6,  # FIXME, not in Table of Optics v3.2
        drawOption=2,
        rTopUpPar=list(mirror),
        rTopDownPar=list(mirror),
        rBotUpPar=list(mirror),
        rBotDownPar=list(mirror),
        rLeftUpPar=list(mirror),
        rLeftDownPar=list(mirror),
        rRightUpPar=list(mirror),
        rRightDownPar=list(mirror),
        rotated=[0, "b_rot", 0],
        relative="arm_guide_4",
        when=when,
    )


def guide_4b(when: str | None = None) -> Guide_gravity:
    """Create guide 4b, the cold guide behind the bender.

    Parameters
    ----------
    when:
        McStas condition for the guide to act, such as ``"bender_in==1"``. None, the default,
        always.

    Returns
    -------
    :
        A new guide.
    """
    guide_4b_in_z = 5.45065  # m, cell B25

    return Guide_gravity(  # row 25
        w1=0.063465631706607162,  # cell E25
        h1=0.044759094882096839,  # cell F25
        w2=GUIDE_4_OUT_WIDTH,
        h2=0.046152195301135549,  # cell H25
        l=0.431,  # cell D25
        mleft=1.5,  # cell O25
        mright=1.5,  # cell P25
        mtop=3.5,  # cell Q25
        mbottom=3.5,  # cell R25
        at=[0, 0, z_along_arm(guide_4b_in_z, GUIDE_4A_IN_Z, GUIDE_4_TILT)],
        relative="arm_guide_4",
        when=when,
    )
