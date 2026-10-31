"""The bispectral extraction guides of T-Rex.

Dimensions are from the Table of Optics V3.2 (the .xlsx file in this folder), whose rows are
named in the comments. The values are the sheet's, not rounded.

The -x wall of the NBOA (Neutron Beam Optics Assembly) is fixed at -47 mm. All of the horizontal
narrowing is taken out of the +x wall, the cold moderator side: the view of the thermal face is
held constant while the cold-side acceptance is trimmed. The guide centre-line is therefore not
the nominal beam axis, but a straight line tilted by half the taper angle.

``Guide_gravity`` is symmetric about its own axis, but two straight walls are reproduced exactly
by a symmetric guide placed on their bisector. Guides 1 to 3 taper at 0.4300 deg (cell N13), so
one arm carries them all and the sections hang off it. Guide 4 tapers at 0.4267 deg (N25) and the
sheet carries the exit width of guide 3 unchanged across the gap, so guide 4 has its own arm,
which keeps the -x wall on the same line. Guides 4a and 4b lie on one axis and share that arm.

Guides 4a (thermal) and 4b (cold, behind the bender) occupy the same slot and share the same exit
plane: they are alternative optics, chosen by the instrument parameter ``bender_in``, not
parallel beamlines. The bender's entrance is oriented toward the cold moderator and its exit
tangent is parallel to the guide axis, so nothing downstream of it is rotated.
"""

import math

from mcdgs.core.components import Arm, Component, Guide_gravity, Pol_bender

BENDER_LENGTH = 0.05  # m, cell B25 minus cell B24
BENDER_RADIUS = 7.2  # m, not in the sheet: from trex_primary.txt
BENDER_ROTATION = -math.degrees(math.atan(BENDER_LENGTH / BENDER_RADIUS))  # −0.39788 deg


def half_taper_angle(w1: float, w2: float, length: float) -> float:
    """Compute the horizontal rotation angle to make the right wall (-x) along z-axis.

    Parameters
    ----------
    w1:
        Width at the entrance, in m.
    w2:
        Width at the exit, in m.
    length:
        Distance between entrance and exit, in m.

    Returns
    -------
    :
        Tilt in deg, negative.
    """
    return -math.degrees(math.atan((w1 - w2) / length / 2))


def entrance_x(width: float, tilt: float) -> float:
    """Compute the x of a guide's entrance centre that puts its -x wall on the fixed wall.

    Parameters
    ----------
    width:
        Width of the guide's entrance, in m.
    tilt:
        Tilt of the guide, from :func:`half_taper_angle`, in deg.

    Returns
    -------
    :
        x of the entrance centre in lab frame, in m.
    """
    nboa_fixed_wall_x = -0.047  # m, the -x wall of the NBOA is fixed at -47 mm
    return width * math.cos(math.radians(tilt)) / 2 + nboa_fixed_wall_x


def z_along_arm(z: float, arm_z: float, tilt: float) -> float:
    """Compute the distance along an arm that places a guide at a lab position.

    Parameters
    ----------
    z:
        z of the guide's entrance in lab frame, in m.
    arm_z:
        z of the arm in lab frame, in m.
    tilt:
        Tilt of the arm, from :func:`half_taper_angle`, in deg.

    Returns
    -------
    :
        Distance in m from the arm, along its tilted z axis.
    """
    return (z - arm_z) / math.cos(math.radians(tilt))


def guides_extraction() -> dict[str, Component]:
    """Create the extraction guides.

    Returns
    -------
    :
        Components by name, in beam order.
    """
    # NBOA constants
    guide_1_in_width = 0.089955  # m, cell E13
    guide_1_in_z = 1.895927  # m, cell B13
    guide_2_in_z = 2.896327  # m, cell B15
    guide_3_in_z = 3.896727  # m, cell B17
    guide_3_out_width = 0.063838  # m, cell G17, the exit of guide 3
    guide_3_out_z = 5.375927  # m, cell C17, the exit of guide 3
    nboa_taper_angle = half_taper_angle(
        guide_1_in_width, guide_3_out_width, guide_3_out_z - guide_1_in_z
    )  # -0.2150 deg, half of cell N13

    # BBGOA constants
    guide_4a_in_z = 5.40065  # m, cell B24, also the entrance of the bender
    guide_4a_in_width = 0.063838  # m, cell E24
    guide_4a_length = 0.481  # m, cell D24
    guide_4b_in_z = 5.45065  # m, cell B25
    guide_4b_in_width = 0.063465631706607162  # m, cell E25
    guide_4b_length = 0.431  # m, cell D25
    guide_4_out_width = 0.060255817017560886  # m, cells G24 and G25
    guide_4_taper_angle = half_taper_angle(
        guide_4a_in_width, guide_4_out_width, guide_4a_length
    )  # -0.2134 deg, half of cells N25 and N26

    # bender constants
    # 63.84mm x 46mm x 50mm (WxHxL), 0.15mm Si channels/ NiTi m=4 curve=7.2m
    bender_mirror = (0.99, 0.0219, 3.02, 4.0, 0.003)  # R0, Qc, alpha, m, W,

    return {
        "arm_NBOA": Arm(
            at=[entrance_x(guide_1_in_width, nboa_taper_angle), 0, guide_1_in_z],
            rotated=[0, nboa_taper_angle, 0],
        ),
        "guide_1": Guide_gravity(  # row 13
            w1=guide_1_in_width,
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
            at=[0, 0, z_along_arm(guide_2_in_z, guide_1_in_z, nboa_taper_angle)],
            relative="arm_NBOA",
        ),
        "guide_3": Guide_gravity(  # row 17
            w1=0.074939,  # cell E17
            h1=0.039164,  # cell F17
            w2=guide_3_out_width,
            h2=0.044511,  # cell H17
            l=1.4792,  # cell D17
            mleft=1.5,  # cell O17
            mright=1.5,  # cell P17
            mtop=4,  # cell Q17, the sheet has 4 to 3.5, marked "issue"
            mbottom=4,  # cell R17
            at=[0, 0, z_along_arm(guide_3_in_z, guide_1_in_z, nboa_taper_angle)],
            relative="arm_NBOA",
        ),
        "arm_guide_4": Arm(
            at=[entrance_x(guide_4a_in_width, guide_4_taper_angle), 0, guide_4a_in_z],
            rotated=[0, guide_4_taper_angle, 0],
        ),
        "guide_4a": Guide_gravity(  # row 24
            w1=guide_4a_in_width,
            h1=0.044593536722227093,  # cell F24
            w2=guide_4_out_width,
            h2=0.046152195301135549,  # cell H24
            l=guide_4a_length,
            mleft=1.5,  # cell O24
            mright=1.5,  # cell P24
            mtop=3.5,  # cell Q24
            mbottom=3.5,  # cell R24
            relative="arm_guide_4",
            when="bender_in==0",
        ),
        "bender": Pol_bender(  # the rest from trex_primary.txt
            xwidth=0.06384,  # not in the sheet: from trex_guide_anyshape.ipynb of mctrex
            yheight=0.046,  # not in the sheet: from trex_guide_anyshape.ipynb of mctrex
            length=BENDER_LENGTH,
            radius=BENDER_RADIUS,
            G=0,  # as every Guide_gravity
            nslit=426,
            d=1e-6,
            drawOption=2,
            rTopUpPar=list(bender_mirror),
            rTopDownPar=list(bender_mirror),
            rBotUpPar=list(bender_mirror),
            rBotDownPar=list(bender_mirror),
            rLeftUpPar=list(bender_mirror),
            rLeftDownPar=list(bender_mirror),
            rRightUpPar=list(bender_mirror),
            rRightDownPar=list(bender_mirror),
            rotated=[0, "b_rot", 0],
            relative="arm_guide_4",
            when="bender_in==1",
        ),
        "guide_4b": Guide_gravity(  # row 25
            w1=guide_4b_in_width,
            h1=0.044759094882096839,  # cell F25
            w2=guide_4_out_width,
            h2=0.046152195301135549,  # cell H25
            l=guide_4b_length,
            mleft=1.5,  # cell O25
            mright=1.5,  # cell P25
            mtop=3.5,  # cell Q25
            mbottom=3.5,  # cell R25
            at=[0, 0, z_along_arm(guide_4b_in_z, guide_4a_in_z, guide_4_taper_angle)],
            relative="arm_guide_4",
            when="bender_in==1",
        ),
    }
