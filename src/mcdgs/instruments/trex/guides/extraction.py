"""The extraction guides 1 to 3 of T-Rex, in the NBOA (Neutron Beam Optics Assembly).

Dimensions are from the Table of Optics V3.2 (the .xlsx file in ``reference/``), whose rows are
named in the comments. The values are the sheet's, not rounded.

The -x wall of the NBOA is fixed at -47 mm. All of the horizontal narrowing is taken out of the +x
wall, the cold moderator side: the view of the thermal face is held constant while the cold-side
acceptance is trimmed. The guide center-line is therefore not the nominal beam axis, but a
straight line tilted by half the taper angle.

``Guide_gravity`` is symmetric about its own axis, but two straight walls are reproduced exactly
by a symmetric guide placed on their bisector. Guides 1 to 3 taper at 0.4300 deg (cell N13), so
one arm carries them all and the sections hang off it. Guide 4 follows on its own arm: see
:mod:`.bispectral`.
"""

import math

from mcdgs.core.components import Arm, Component, Guide_gravity

# the entrance of guide 1, row 13: the source focuses on it, and the source monitor is it
GUIDE_1_IN_Z = 1.895927  # m, cell B13
GUIDE_1_IN_WIDTH = 0.089955  # m, cell E13
GUIDE_1_IN_HEIGHT = 0.029992  # m, cell F13

# x of the fixed wall of the NBOA, which the -x wall of the guides is on: the source focuses
# up to it on both sides
NBOA_FIXED_WALL_X = -0.047  # FIXME, not in Table of Optics v3.2, in m


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
    """Compute the x of a guide's entrance center that puts its -x wall on the fixed wall.

    Parameters
    ----------
    width:
        Width of the guide's entrance, in m.
    tilt:
        Tilt of the guide, from :func:`half_taper_angle`, in deg.

    Returns
    -------
    :
        x of the entrance center in lab frame, in m.
    """
    return width * math.cos(math.radians(tilt)) / 2 + NBOA_FIXED_WALL_X


def z_along_arm(z_lab: float, arm_z_lab: float, tilt: float) -> float:
    """Compute the distance along an arm that places a guide at a lab position.

    Parameters
    ----------
    z_lab:
        z of the guide's entrance in lab frame, in m.
    arm_z_lab:
        z of the arm in lab frame, in m.
    tilt:
        Tilt of the arm, from :func:`half_taper_angle`, in deg.

    Returns
    -------
    :
        Distance in m from the arm, along its tilted z axis.
    """
    return (z_lab - arm_z_lab) / math.cos(math.radians(tilt))


def guides_extraction() -> dict[str, Component]:
    """Create the extraction guides 1 to 3. Guide 4 follows them: see :mod:`.bispectral`.

    Returns
    -------
    :
        Components by name, in beam order.
    """
    # NBOA constants
    guide_2_in_z = 2.896327  # m, cell B15
    guide_3_in_z = 3.896727  # m, cell B17
    guide_3_out_width = 0.063838  # m, cell G17, the exit of guide 3
    guide_3_out_z = 5.375927  # m, cell C17, the exit of guide 3
    nboa_taper_angle = half_taper_angle(
        GUIDE_1_IN_WIDTH, guide_3_out_width, guide_3_out_z - GUIDE_1_IN_Z
    )  # -0.2150 deg, half of cell N13

    return {
        "arm_nboa": Arm(
            at=[entrance_x(GUIDE_1_IN_WIDTH, nboa_taper_angle), 0, GUIDE_1_IN_Z],
            rotated=[0, nboa_taper_angle, 0],
        ),
        "guide_1": Guide_gravity(  # row 13
            w1=GUIDE_1_IN_WIDTH,
            h1=GUIDE_1_IN_HEIGHT,
            w2=0.082455,  # cell G13
            h2=0.034937,  # cell H13
            l=0.9994,  # cell D13
            mleft=1.5,  # cell O13
            mright=1.5,  # cell P13
            mtop=4,  # cell Q13
            mbottom=4,  # cell R13
            relative="arm_nboa",
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
            at=[0, 0, z_along_arm(guide_2_in_z, GUIDE_1_IN_Z, nboa_taper_angle)],
            relative="arm_nboa",
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
            at=[0, 0, z_along_arm(guide_3_in_z, GUIDE_1_IN_Z, nboa_taper_angle)],
            relative="arm_nboa",
        ),
    }
