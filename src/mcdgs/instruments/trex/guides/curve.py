"""The curved guides of T-Rex, from guide 5 to guide 29.

Dimensions are from the Table of Optics V3.2 (the .xlsx file in ``reference/``), whose rows are
named in the comments. The values are the sheet's, not rounded. The m-values are those at the
entrance, which are the same as at the exit.

The sheet's positions are distances along the beam axis. The distance between two entrances is
therefore the difference of their positions, while lab z falls behind the sheet where the axis
turns.

The curve, 12 km in radius, starts at guide 5 and ends at the exit of guide 29. It is made of
straight sections, each at the angle the arc has at the section's midpoint, so a section turns
from the one before it by the arc between their midpoints, gaps included.

``arm_start_of_curve``, at the entrance of guide 5, takes out the tilt of guide 4: the fixed -x
wall ends there, and from guide 5 on each guide is symmetric about the beam axis. Guide 5 turns
from it by half its own length. Guide 29 lies at the angle the arc has at its midpoint, so
``arm_end_of_curve``, at its exit, takes up the arc over the second half of its length, for the
axis to leave the curve on its tangent.

The gaps for the choppers and the beam monitor are arms in the chain, for what goes in them to be
placed relative to.
"""

from mcdgs.core.components import Arm, Component

from .bispectral import GUIDE_4_TILT, GUIDE_4A_IN_Z
from .extraction import z_along_arm
from .sections import arc_angle, chain, gap, guide, split_after

CURVE_RADIUS = 12000.0  # m, cell N32
GUIDE_29_OUT_Z = 95.99955  # m, cell C106
ARM_BWC1_Z = 31.94904  # m, cell B65
BWC1_GAP_LENGTH = 0.03  # m, cell D65
ARM_BWC2_Z = 39.94004  # m, cell B73
BWC2_GAP_LENGTH = 0.03  # m, cell D73


def turn(distance: float) -> float:
    """Compute the turn of the beam axis along a distance of the curve.

    Parameters
    ----------
    distance:
        Distance along the curve, in m.

    Returns
    -------
    :
        Rotation about +y, in deg.
    """
    return arc_angle(distance, CURVE_RADIUS)


def guides_curve() -> list[dict[str, Component]]:
    """Create the curved guides.

    Returns
    -------
    :
        Components by name, in beam order, in three pieces: up to ``arm_bwc1``, up to
        ``arm_bwc2``, and the rest.
    """
    guide_5_in_z = 5.916  # m, cell B32
    guide_5_length = 2.7625  # m, cell D32
    guide_29_length = 1.998  # m, cell D106

    sections = {
        "guide_5": guide(  # row 32
            z=guide_5_in_z,
            w1=0.06,  # cell E32
            h1=0.046260695185495173,  # cell F32
            w2=0.06,  # cell G32
            h2=0.053982339063629709,  # cell H32
            l=guide_5_length,
            mleft=1.5,  # cell O32
            mright=1.5,  # cell P32
            mtop=3.5,  # cell Q32
            mbottom=3.5,  # cell R32
            rotated=[0, turn(guide_5_length / 2), 0],
            relative="arm_start_of_curve",
        ),
        "guide_6": guide(  # row 38
            z=8.6915,  # cell B38
            w1=0.06,  # cell E38
            h1=0.05401466907725186,  # cell F38
            w2=0.06,  # cell G38
            h2=0.060228315434988687,  # cell H38
            l=2.7625,  # cell D38
            mleft=1.5,  # cell O38
            mright=1.5,  # cell P38
            mtop=3,  # cell Q38
            mbottom=3,  # cell R38
        ),
        "guide_7": guide(  # row 44
            z=11.467,  # cell B44
            w1=0.06,  # cell E44
            h1=0.060254788378139509,  # cell F44
            w2=0.06,  # cell G44
            h2=0.067415573123614891,  # cell H44
            l=4,  # cell D44
            mleft=1.5,  # cell O44
            mright=1.5,  # cell P44
            mtop=3,  # cell Q44
            mbottom=3,  # cell R44
        ),
        "guide_8": guide(  # row 46
            z=15.469,  # cell B46
            w1=0.06,  # cell E46
            h1=0.067418714855804978,  # cell F46
            w2=0.06,  # cell G46
            h2=0.072976133271138181,  # cell H46
            l=4,  # cell D46
            mleft=1.5,  # cell O46
            mright=1.5,  # cell P46
            mtop=3,  # cell Q46
            mbottom=3,  # cell R46
        ),
        "guide_9": guide(  # row 48
            z=19.471,  # cell B48
            w1=0.06,  # cell E48
            h1=0.072978576629709963,  # cell F48
            w2=0.06,  # cell G48
            h2=0.0765252384,  # cell H48
            l=3.2167,  # cell D48
            mleft=1.5,  # cell O48
            mright=1.5,  # cell P48
            mtop=3,  # cell Q48
            mbottom=3,  # cell R48
        ),
        "guide_10": guide(  # row 54
            z=22.7007,  # cell B54
            w1=0.06,  # cell E54
            h1=0.0765380906,  # cell F54
            w2=0.06,  # cell G54
            h2=0.078023975583694863,  # cell H54
            l=1.588,  # cell D54
            mleft=2,  # cell O54
            mright=2,  # cell P54
            mtop=2.5,  # cell Q54
            mbottom=2.5,  # cell R54
        ),
        "guide_11": guide(  # row 60
            z=24.3017,  # cell B60
            w1=0.06,  # cell E60
            h1=0.078035464780223663,  # cell F60
            w2=0.06,  # cell G60
            h2=0.080841682168482322,  # cell H60
            l=3.6349,  # cell D60
            mleft=1.5,  # cell O60
            mright=2.5,  # cell P60
            mtop=2.5,  # cell Q60
            mbottom=2.5,  # cell R60
        ),
        "guide_12": guide(  # row 62
            z=27.94614,  # cell B62
            w1=0.06,  # cell E62
            h1=0.080848019233867547,  # cell F62
            w2=0.06,  # cell G62
            h2=0.082063226084044189,  # cell H62
            l=1.9998,  # cell D62
            mleft=1.5,  # cell O62
            mright=2.5,  # cell P62
            mtop=2.5,  # cell Q62
            mbottom=2.5,  # cell R62
        ),
        "guide_13a": guide(  # row 64
            z=29.94924,  # cell B64
            w1=0.06,  # cell E64
            h1=0.082065047956342354,  # cell F64
            w2=0.06,  # cell G64
            h2=0.08306090771746841,  # cell H64
            l=1.9998,  # cell D64
            mleft=1.5,  # cell O64
            mright=2.5,  # cell P64
            mtop=2.5,  # cell Q64
            mbottom=2.5,  # cell R64
        ),
        "arm_bwc1": gap(  # row 65, bandwidth chopper 1
            z=ARM_BWC1_Z,
            length=BWC1_GAP_LENGTH,
        ),
        "guide_13b": guide(  # row 66
            z=31.97904,  # cell B66
            w1=0.06,  # cell E66
            h1=0.083074222418113493,  # cell F66
            w2=0.06,  # cell G66
            h2=0.083855788380960263,  # cell H66
            l=1.9998,  # cell D66
            mleft=1.5,  # cell O66
            mright=2.5,  # cell P66
            mtop=2.5,  # cell Q66
            mbottom=2.5,  # cell R66
        ),
        "guide_14a": guide(  # row 68
            z=33.98214,  # cell B68
            w1=0.06,  # cell E68
            h1=0.083856907127152837,  # cell F68
            w2=0.06,  # cell G68
            h2=0.084432926758712085,  # cell H68
            l=1.9998,  # cell D68
            mleft=1.5,  # cell O68
            mright=2.5,  # cell P68
            mtop=2.5,  # cell Q68
            mbottom=2.5,  # cell R68
        ),
        "guide_14b": guide(  # row 70
            z=35.98264,  # cell B70
            w1=0.06,  # cell E70
            h1=0.084433093,  # cell F70
            w2=0.06,  # cell G70
            h2=0.084802960774443392,  # cell H70
            l=1.9633,  # cell D70
            mleft=1.5,  # cell O70
            mright=2.5,  # cell P70
            mtop=2.5,  # cell Q70
            mbottom=2.5,  # cell R70
        ),
        "guide_15a": guide(  # row 72
            z=37.94924,  # cell B72
            w1=0.06,  # cell E72
            h1=0.084803420972849011,  # cell F72
            w2=0.06,  # cell G72
            h2=0.084982900353593493,  # cell H72
            l=1.9908,  # cell D72
            mleft=1.5,  # cell O72
            mright=2.5,  # cell P72
            mtop=2.5,  # cell Q72
            mbottom=2.5,  # cell R72
        ),
        "arm_bwc2": gap(  # row 73, bandwidth chopper 2
            z=ARM_BWC2_Z,
            length=BWC2_GAP_LENGTH,
        ),
        "guide_15b": guide(  # row 74
            z=39.97004,  # cell B74
            w1=0.06,  # cell E74
            h1=0.0849841094,  # cell F74
            w2=0.06,  # cell G74
            h2=0.0849646725,  # cell H74
            l=2.0012,  # cell D74
            mleft=1.5,  # cell O74
            mright=2.5,  # cell P74
            mtop=2.5,  # cell Q74
            mbottom=2.5,  # cell R74
        ),
        "arm_bm1": gap(  # row 77, beam monitor 1
            z=41.97786,  # cell B77
            length=0.02,  # cell D77
        ),
        "guide_16": guide(  # row 80
            z=42.00155,  # cell B80
            w1=0.06,  # cell E80
            h1=0.085,  # cell F80
            w2=0.06,  # cell G80
            h2=0.085,  # cell H80
            l=3.9967,  # cell D80
            mleft=1.5,  # cell O80
            mright=2.5,  # cell P80
            mtop=2.5,  # cell Q80
            mbottom=2.5,  # cell R80
        ),
        "guide_17": guide(  # row 82
            z=46.00155,  # cell B82
            w1=0.06,  # cell E82
            h1=0.085,  # cell F82
            w2=0.06,  # cell G82
            h2=0.085,  # cell H82
            l=3.9967,  # cell D82
            mleft=1.5,  # cell O82
            mright=2.5,  # cell P82
            mtop=2.5,  # cell Q82
            mbottom=2.5,  # cell R82
        ),
        "guide_18": guide(  # row 84
            z=50.00155,  # cell B84
            w1=0.06,  # cell E84
            h1=0.085,  # cell F84
            w2=0.06,  # cell G84
            h2=0.085,  # cell H84
            l=3.9967,  # cell D84
            mleft=1.5,  # cell O84
            mright=2.5,  # cell P84
            mtop=2.5,  # cell Q84
            mbottom=2.5,  # cell R84
        ),
        "guide_19": guide(  # row 86
            z=54.00155,  # cell B86
            w1=0.06,  # cell E86
            h1=0.085,  # cell F86
            w2=0.06,  # cell G86
            h2=0.085,  # cell H86
            l=3.9967,  # cell D86
            mleft=1.5,  # cell O86
            mright=2.5,  # cell P86
            mtop=2.5,  # cell Q86
            mbottom=2.5,  # cell R86
        ),
        "guide_20": guide(  # row 88
            z=58.00155,  # cell B88
            w1=0.06,  # cell E88
            h1=0.085,  # cell F88
            w2=0.06,  # cell G88
            h2=0.085,  # cell H88
            l=3.9967,  # cell D88
            mleft=1.5,  # cell O88
            mright=2.5,  # cell P88
            mtop=2.5,  # cell Q88
            mbottom=2.5,  # cell R88
        ),
        "guide_21": guide(  # row 90
            z=62.00155,  # cell B90
            w1=0.06,  # cell E90
            h1=0.085,  # cell F90
            w2=0.06,  # cell G90
            h2=0.085,  # cell H90
            l=3.9967,  # cell D90
            mleft=1.5,  # cell O90
            mright=2.5,  # cell P90
            mtop=2.5,  # cell Q90
            mbottom=2.5,  # cell R90
        ),
        "guide_22": guide(  # row 92
            z=66.00155,  # cell B92
            w1=0.06,  # cell E92
            h1=0.085,  # cell F92
            w2=0.06,  # cell G92
            h2=0.085,  # cell H92
            l=3.9967,  # cell D92
            mleft=1.5,  # cell O92
            mright=2.5,  # cell P92
            mtop=2.5,  # cell Q92
            mbottom=2.5,  # cell R92
        ),
        "guide_23": guide(  # row 94
            z=70.00155,  # cell B94
            w1=0.06,  # cell E94
            h1=0.085,  # cell F94
            w2=0.06,  # cell G94
            h2=0.085,  # cell H94
            l=3.9967,  # cell D94
            mleft=1.5,  # cell O94
            mright=2.5,  # cell P94
            mtop=2.5,  # cell Q94
            mbottom=2.5,  # cell R94
        ),
        "guide_24": guide(  # row 96
            z=74.00155,  # cell B96
            w1=0.06,  # cell E96
            h1=0.085,  # cell F96
            w2=0.06,  # cell G96
            h2=0.085,  # cell H96
            l=3.9967,  # cell D96
            mleft=1.5,  # cell O96
            mright=2.5,  # cell P96
            mtop=2.5,  # cell Q96
            mbottom=2.5,  # cell R96
        ),
        "guide_25": guide(  # row 98
            z=78.00155,  # cell B98
            w1=0.06,  # cell E98
            h1=0.085,  # cell F98
            w2=0.06,  # cell G98
            h2=0.085,  # cell H98
            l=3.9967,  # cell D98
            mleft=1.5,  # cell O98
            mright=2.5,  # cell P98
            mtop=2.5,  # cell Q98
            mbottom=2.5,  # cell R98
        ),
        "guide_26": guide(  # row 100
            z=82.00155,  # cell B100
            w1=0.06,  # cell E100
            h1=0.085,  # cell F100
            w2=0.06,  # cell G100
            h2=0.085,  # cell H100
            l=3.9967,  # cell D100
            mleft=1,  # cell O100
            mright=2.5,  # cell P100
            mtop=2.5,  # cell Q100
            mbottom=2.5,  # cell R100
        ),
        "guide_27": guide(  # row 102
            z=86.00155,  # cell B102
            w1=0.06,  # cell E102
            h1=0.085,  # cell F102
            w2=0.06,  # cell G102
            h2=0.085,  # cell H102
            l=3.9967,  # cell D102
            mleft=1,  # cell O102
            mright=2.5,  # cell P102
            mtop=2.5,  # cell Q102
            mbottom=2.5,  # cell R102
        ),
        "guide_28": guide(  # row 104
            z=90.00155,  # cell B104
            w1=0.06,  # cell E104
            h1=0.085,  # cell F104
            w2=0.06,  # cell G104
            h2=0.085,  # cell H104
            l=3.9967,  # cell D104
            mleft=1,  # cell O104
            mright=2.5,  # cell P104
            mtop=2.5,  # cell Q104
            mbottom=2.5,  # cell R104
        ),
        "guide_29": guide(  # row 106
            z=94.00155,  # cell B106
            w1=0.06,  # cell E106
            h1=0.085,  # cell F106
            w2=0.06,  # cell G106
            h2=0.085,  # cell H106
            l=guide_29_length,
            mleft=1,  # cell O106
            mright=2.5,  # cell P106
            mtop=2.5,  # cell Q106
            mbottom=2.5,  # cell R106
        ),
    }
    components = {
        "arm_start_of_curve": Arm(
            at=[0, 0, z_along_arm(guide_5_in_z, GUIDE_4A_IN_Z, GUIDE_4_TILT)],
            rotated=[0, -GUIDE_4_TILT, 0],
            relative="arm_guide_4",
        ),
        **chain(sections, CURVE_RADIUS),
        "arm_end_of_curve": Arm(
            at=[0, 0, guide_29_length],
            rotated=[0, turn(guide_29_length / 2), 0],
            relative="guide_29",
        ),
    }
    return split_after(components, "arm_bwc1", "arm_bwc2")
