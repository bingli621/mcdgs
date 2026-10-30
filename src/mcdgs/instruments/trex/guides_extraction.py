"""The bispectral extraction guides of T-Rex.

Dimensions are from the Table of Optics V3.2 (the .xlsx file in this folder), whose rows are
named in the comments. The values are the sheet's, not rounded.

The -x wall of the NBOA (Neutron Beam Optics Assembly) is fixed at -47 mm. All of the horizontal
narrowing is taken out of the +x wall, the cold moderator side: the view of the thermal face is
held constant while the cold-side acceptance is trimmed. The guide centre-line is therefore not
the nominal beam axis, but a straight line tilted by half the taper angle.

``Guide_gravity`` is symmetric about its own axis, but two straight walls are reproduced exactly
by a symmetric guide placed on their bisector. The tilt is taken as the same for every section,
so one arm carries the whole NBOA and the sections hang off it. The tilt is that of guides 1 to
3, which taper at 0.4300 deg (cell N13). This is an approximation for guide 4, which tapers at
0.4267 deg (N25).

Guides 4a (thermal) and 4b (cold, behind the bender) occupy the same slot and share the same exit
plane: they are alternative optics, chosen by the instrument parameter ``bender``, not parallel
beamlines. The bender's entrance is oriented toward the cold moderator and its exit tangent lies
on the guide axis, so nothing downstream of it is rotated.
"""

import math

from mcdgs.core.components import Arm, Component, Guide_gravity, Pol_bender

GUIDE_1_IN_WIDTH = 0.089955  # m, cell E13
GUIDE_1_IN_Z = 1.895927  # m, cell B13
GUIDE_3_OUT_WIDTH = 0.063838  # m, cell G17, the exit of guide 3
GUIDE_3_OUT_Z = 5.375927  # m, cell C17, the exit of guide 3
NBOA_FIXED_WALL_X = -0.047  # m, not in the sheet: from guides_extraction.py of mctrex

GUIDE_2_IN_Z = 2.896327  # m, cell B15
GUIDE_3_IN_Z = 3.896727  # m, cell B17
GUIDE_4A_IN_Z = 5.40065  # m, cell B24, also the entrance of the bender
GUIDE_4B_IN_Z = 5.45065  # m, cell B25


NBOA_TAPER_ANGLE = math.degrees(
    math.atan((GUIDE_1_IN_WIDTH - GUIDE_3_OUT_WIDTH) / (GUIDE_3_OUT_Z - GUIDE_1_IN_Z) / 2)
)  # 0.2150 deg, half of cell N13

NBOA_IN_X = GUIDE_1_IN_WIDTH * math.cos(math.radians(NBOA_TAPER_ANGLE)) / 2 + NBOA_FIXED_WALL_X

BENDER_LENGTH = 0.05  # m, cell B25 minus cell B24
BENDER_RADIUS = 7.2  # m, not in the sheet: from trex_primary.txt
BENDER_ROTATION = -math.degrees(math.atan(BENDER_LENGTH / BENDER_RADIUS))  # −0.39788 deg
BENDER_MIRROR = (0.99, 0.0219, 3.02, 4.0, 0.003)  # R0, Qc, alpha, m, W, from trex_primary.txt


def z_along_arm(z: float) -> float:
    """Compute the distance along ``arm_NBOA`` that places a guide at a lab position.

    Parameters
    ----------
    z:
        z of the guide's entrance in lab frame, in m.

    Returns
    -------
    :
        Distance in m from the arm, along its tilted z axis.
    """
    return (z - GUIDE_1_IN_Z) / math.cos(math.radians(NBOA_TAPER_ANGLE))


def guides_extraction() -> dict[str, Component]:
    """Create the extraction guides.

    Returns
    -------
    :
        New components by name, in beam order.
    """
    return {
        "arm_NBOA": Arm(at=[NBOA_IN_X, 0, GUIDE_1_IN_Z], rotated=[0, -NBOA_TAPER_ANGLE, 0]),
        "guide_1": Guide_gravity(  # row 13
            w1=GUIDE_1_IN_WIDTH,
            h1=0.029992,  # cell F13
            w2=0.082455,  # cell G13
            h2=0.034937,  # cell H13
            l=0.9994,  # cell D13
            mleft=1.5,  # cell O13
            mright=1.5,  # cell P13
            mtop=4,  # cell Q13
            mbottom=4,  # cell R13
            relative="arm_NBOA",
        ),
        "guide_2": Guide_gravity(  # row 15
            w1=0.082447,  # cell E15
            h1=0.034941,  # cell F15
            w2=0.074948,  # cell G15
            h2=0.03916,  # cell H15
            l=0.9992,  # cell D15
            mleft=1.5,  # cell O15
            mright=1.5,  # cell P15
            mtop=4,  # cell Q15
            mbottom=4,  # cell R15
            at=[0, 0, z_along_arm(GUIDE_2_IN_Z)],
            relative="arm_NBOA",
        ),
        "guide_3": Guide_gravity(  # row 17
            w1=0.074939,  # cell E17
            h1=0.039164,  # cell F17
            w2=GUIDE_3_OUT_WIDTH,
            h2=0.044511,  # cell H17
            l=1.4792,  # cell D17
            mleft=1.5,  # cell O17
            mright=1.5,  # cell P17
            mtop=4,  # cell Q17, the sheet has 4 to 3.5, marked "issue"
            mbottom=4,  # cell R17
            at=[0, 0, z_along_arm(GUIDE_3_IN_Z)],
            relative="arm_NBOA",
        ),
        "guide_4a": Guide_gravity(  # row 24
            w1=0.063838,  # cell E24
            h1=0.044593536722227093,  # cell F24
            w2=0.060255817017560886,  # cell G24
            h2=0.046152195301135549,  # cell H24
            l=0.481,  # cell D24
            mleft=1.5,  # cell O24
            mright=1.5,  # cell P24
            mtop=3.5,  # cell Q24
            mbottom=3.5,  # cell R24
            at=[0, 0, z_along_arm(GUIDE_4A_IN_Z)],
            relative="arm_NBOA",
            when="bender==0",
        ),
        "bender": Pol_bender(  # size of row 25, the rest from trex_primary.txt
            xwidth=0.063466,  # estimated using cell E25
            yheight=0.044759,  # estimated using cell F25
            length=BENDER_LENGTH,
            radius=BENDER_RADIUS,
            G=0,  # as every Guide_gravity
            nslit=426,
            d=1e-6,
            drawOption=2,
            rTopUpPar=list(BENDER_MIRROR),
            rTopDownPar=list(BENDER_MIRROR),
            rBotUpPar=list(BENDER_MIRROR),
            rBotDownPar=list(BENDER_MIRROR),
            rLeftUpPar=list(BENDER_MIRROR),
            rLeftDownPar=list(BENDER_MIRROR),
            rRightUpPar=list(BENDER_MIRROR),
            rRightDownPar=list(BENDER_MIRROR),
            at=[0, 0, z_along_arm(GUIDE_4A_IN_Z)],
            rotated=[0, "b_rot", 0],
            relative="arm_NBOA",
            when="bender==1",
        ),
        "guide_4b": Guide_gravity(  # row 25
            w1=0.063465631706607162,  # cell E25
            h1=0.044759094882096839,  # cell F25
            w2=0.060255817017560886,  # cell G25
            h2=0.046152195301135549,  # cell H25
            l=0.431,  # cell D25
            mleft=1.5,  # cell O25
            mright=1.5,  # cell P25
            mtop=3.5,  # cell Q25
            mbottom=3.5,  # cell R25
            at=[0, 0, z_along_arm(GUIDE_4B_IN_Z)],
            relative="arm_NBOA",
            when="bender==1",
        ),
    }
