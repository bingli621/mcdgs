"""The straight guides of CSPEC, rows 78 to 189, and the sample position.

Dimensions are from the Specification sheet of CSPEC nl-specs 1.E (the .xlsx file in
``reference/``), whose rows and cells are named in the comments, in mm and deg as the sheet gives
them. The values are the sheet's, not rounded. Each guide and gap is placed by its own start in
ISCS: see :mod:`.frame`. The gaps of 1 to 10 mm between the elements are left empty.

The trajectory is straight from row 77 on, 422 mm below the beamline and turned 0.745 deg from
the insert, so every element has the same angles. The sheet approximates the ellipses of the
guide with straight tapered elements: diverging from row 78, converging in height from row 117,
and converging in width from row 139 and in height from row 153.

The guide exchanger of row 191 is not modeled: the sheet leaves it to be defined. The beam flies
free from the exit of W03-09-02 to the sample position.

The gaps of bandwidth chopper 3, the pulse-shaping choppers and the monochromatic choppers are
arms, for the choppers to be placed relative to.
"""

import math

from mcdgs.core.components import Arm, Component

from .frame import MM_PER_M, arm, guide

# the flight path from the moderator to the start of each chopper gap, along the trajectory
ARM_BWC3_PATH = 104.59713344998615  # m, cell U129
ARM_PSC_PATH = 105.65601344998615  # m, cell U133
ARM_MC_PATH = 158.45582344998624  # m, cell U187

# the exit of W03-09-02, the last guide
GUIDE_09_02_OUT_PATH = 158.82879344998625  # m, cell U189
GUIDE_09_02_OUT_X = 158791.60693599467  # mm, cell K189
GUIDE_09_02_OUT_Y = 1695.5697820437094  # mm, cell L189
GUIDE_09_02_OUT_Z = -422.42085644514486  # mm, cell M189
GUIDE_09_02_PHI = -0.7448451336700711  # deg, cell G189

SAMPLE_X = 159966.76  # mm, cell B35 of Definitions
SAMPLE_Y = 1710.64  # mm, cell C35 of Definitions
SAMPLE_Z = -422.42  # mm, cell D35 of Definitions

# the flight path from the moderator to the sample: along the trajectory to the exit of the last
# guide, then straight on
SAMPLE_PATH = (
    GUIDE_09_02_OUT_PATH
    + math.dist(
        (GUIDE_09_02_OUT_X, GUIDE_09_02_OUT_Y, GUIDE_09_02_OUT_Z),
        (SAMPLE_X, SAMPLE_Y, SAMPLE_Z),
    )
    / MM_PER_M
)  # 160.004 m


def guides_straight() -> list[dict[str, Component]]:
    """Create the straight guides.

    Returns
    -------
    :
        Components by name, in beam order, in four pieces: up to ``arm_bwc3``, up to
        ``arm_psc``, up to ``arm_mc``, and the rest.
    """
    to_bwc3: dict[str, Component] = {
        "guide_06_14": guide(  # row 78, diverging in width and height, to row 90
            x=54376.69952720585,  # cell H78
            y=338.0995282097471,  # cell I78
            z=-422.42085644514486,  # cell J78
            chi=0,  # cell F78
            phi=-0.7448451336700711,  # cell G78
            length=2000,  # cell T78
            width_in=69.99990607313028,  # cell X78
            width_out=84.10595610236263,  # cell Y78
            height_in=55.002905219987134,  # cell Z78
            height_out=68.36992967200929,  # cell AA78
            mtop=2,  # cell AB78
            mbottom=2,  # cell AC78
            mleft=2,  # cell AD78
            mright=2,  # cell AE78
        ),
        "guide_06_15": guide(  # row 80
            x=56378.5303605883,  # cell H80
            y=364.12479515027485,  # cell I80
            z=-422.42085644514486,  # cell J80
            chi=0,  # cell F80
            phi=-0.7448451336700711,  # cell G80
            length=2400,  # cell T80
            width_in=84.11809754957477,  # cell X80
            width_out=96.63195215328302,  # cell Y80
            height_in=68.38160722903005,  # cell Z80
            height_out=80.72912641027675,  # cell AA80
            mtop=2,  # cell AB80
            mbottom=2,  # cell AC80
            mleft=2,  # cell AD80
            mright=2,  # cell AE80
        ),
        "guide_06_16": guide(  # row 82
            x=58779.327478945575,  # cell H82
            y=395.3369159915372,  # cell I82
            z=-422.42085644514486,  # cell J82
            chi=0,  # cell F82
            phi=-0.7448451336700711,  # cell G82
            length=2400,  # cell T82
            width_in=96.6364210217418,  # cell X82
            width_out=105.93599257907043,  # cell Y82
            height_in=80.7336837495513,  # cell Z82
            height_out=90.60251606551464,  # cell AA82
            mtop=2,  # cell AB82
            mbottom=2,  # cell AC82
            mleft=2,  # cell AD82
            mright=2,  # cell AE82
        ),
        "guide_06_17": guide(  # row 84
            x=61181.12451280404,  # cell H84
            y=426.5620364666359,  # cell I84
            z=-422.42085644514486,  # cell J84
            chi=0,  # cell F84
            phi=-0.7448451336700711,  # cell G84
            length=1600,  # cell T84
            width_in=105.94265809386454,  # cell X84
            width_out=110.76818790974775,  # cell Y84
            height_in=90.6099462352587,  # cell Z84
            height_out=96.19716754892261,  # cell AA84
            mtop=2,  # cell AB84
            mbottom=2,  # cell AC84
            mleft=2,  # cell AD84
            mright=2,  # cell AE84
        ),
        "guide_06_18": guide(  # row 86
            x=62781.989230209285,  # cell H86
            y=447.37445023875625,  # cell I86
            z=-422.42085644514486,  # cell J86
            chi=0,  # cell F86
            phi=-0.7448451336700711,  # cell G86
            length=1600,  # cell T86
            width_in=110.77090108093587,  # cell X86
            width_out=114.66421844673911,  # cell Y86
            height_in=96.2004496021494,  # cell Z86
            height_out=101.15069609555333,  # cell AA86
            mtop=2,  # cell AB86
            mbottom=2,  # cell AC86
            mleft=2,  # cell AD86
            mright=2,  # cell AE86
        ),
        "guide_06_19": guide(  # row 88
            x=64383.85386311572,  # cell H88
            y=468.199863644713,  # cell I88
            z=-422.42085644514486,  # cell J88
            chi=0,  # cell F88
            phi=-0.7448451336700711,  # cell G88
            length=2000,  # cell T88
            width_in=114.66,  # cell X88
            width_out=114.66,  # cell Y88
            height_in=101.15652583072013,  # cell Z88
            height_out=106.58253334116885,  # cell AA88
            mtop=2,  # cell AB88
            mbottom=2,  # cell AC88
            mleft=2,  # cell AD88
            mright=2,  # cell AE88
        ),
        "guide_06_20": guide(  # row 90
            x=66384.68478099698,  # cell H90
            y=494.21213095140433,  # cell I90
            z=-422.42085644514486,  # cell J90
            chi=0,  # cell F90
            phi=-0.7448451336700711,  # cell G90
            length=2000,  # cell T90
            width_in=114.66,  # cell X90
            width_out=114.66,  # cell Y90
            height_in=106.58505486882994,  # cell Z90
            height_out=111.27988185917508,  # cell AA90
            mtop=2,  # cell AB90
            mbottom=2,  # cell AC90
            mleft=2,  # cell AD90
            mright=2,  # cell AE90
        ),
        "guide_06_21": guide(  # row 92, constant, to row 116
            x=68386.51561437943,  # cell H92
            y=520.2373978919321,  # cell I92
            z=-422.42085644514486,  # cell J92
            chi=0,  # cell F92
            phi=-0.7448451336700711,  # cell G92
            length=2000,  # cell T92
            width_in=114.66,  # cell X92
            width_out=114.66,  # cell Y92
            height_in=111.27988185917508,  # cell Z92
            height_out=111.27988185917508,  # cell AA92
            mtop=2,  # cell AB92
            mbottom=2,  # cell AC92
            mleft=2,  # cell AD92
            mright=2,  # cell AE92
        ),
        "guide_06_22": guide(  # row 94
            x=70387.3465322607,  # cell H94
            y=546.2496651986235,  # cell I94
            z=-422.42085644514486,  # cell J94
            chi=0,  # cell F94
            phi=-0.7448451336700711,  # cell G94
            length=2000,  # cell T94
            width_in=114.66,  # cell X94
            width_out=114.66,  # cell Y94
            height_in=111.27988185917508,  # cell Z94
            height_out=111.27988185917508,  # cell AA94
            mtop=2,  # cell AB94
            mbottom=2,  # cell AC94
            mleft=2,  # cell AD94
            mright=2,  # cell AE94
        ),
        "guide_06_23": guide(  # row 96
            x=72389.17736564315,  # cell H96
            y=572.2749321391512,  # cell I96
            z=-422.42085644514486,  # cell J96
            chi=0,  # cell F96
            phi=-0.7448451336700711,  # cell G96
            length=2000,  # cell T96
            width_in=114.66,  # cell X96
            width_out=114.66,  # cell Y96
            height_in=111.27988185917508,  # cell Z96
            height_out=111.27988185917508,  # cell AA96
            mtop=2,  # cell AB96
            mbottom=2,  # cell AC96
            mleft=2,  # cell AD96
            mright=2,  # cell AE96
        ),
        "guide_06_24": guide(  # row 98
            x=74390.00828352441,  # cell H98
            y=598.2871994458426,  # cell I98
            z=-422.42085644514486,  # cell J98
            chi=0,  # cell F98
            phi=-0.7448451336700711,  # cell G98
            length=2000,  # cell T98
            width_in=114.66,  # cell X98
            width_out=114.66,  # cell Y98
            height_in=111.27988185917508,  # cell Z98
            height_out=111.27988185917508,  # cell AA98
            mtop=2,  # cell AB98
            mbottom=2,  # cell AC98
            mleft=2,  # cell AD98
            mright=2,  # cell AE98
        ),
        "guide_06_25": guide(  # row 100
            x=76391.83911690686,  # cell H100
            y=624.3124663863704,  # cell I100
            z=-422.42085644514486,  # cell J100
            chi=0,  # cell F100
            phi=-0.7448451336700711,  # cell G100
            length=2000,  # cell T100
            width_in=114.66,  # cell X100
            width_out=114.66,  # cell Y100
            height_in=111.27988185917508,  # cell Z100
            height_out=111.27988185917508,  # cell AA100
            mtop=2,  # cell AB100
            mbottom=2,  # cell AC100
            mleft=2,  # cell AD100
            mright=2,  # cell AE100
        ),
        "guide_06_26": guide(  # row 102
            x=78392.67003478813,  # cell H102
            y=650.3247336930617,  # cell I102
            z=-422.42085644514486,  # cell J102
            chi=0,  # cell F102
            phi=-0.7448451336700711,  # cell G102
            length=2000,  # cell T102
            width_in=114.66,  # cell X102
            width_out=114.66,  # cell Y102
            height_in=111.27988185917508,  # cell Z102
            height_out=111.27988185917508,  # cell AA102
            mtop=2,  # cell AB102
            mbottom=2,  # cell AC102
            mleft=2,  # cell AD102
            mright=2,  # cell AE102
        ),
        "guide_06_27": guide(  # row 104
            x=80394.50086817058,  # cell H104
            y=676.3500006335895,  # cell I104
            z=-422.42085644514486,  # cell J104
            chi=0,  # cell F104
            phi=-0.7448451336700711,  # cell G104
            length=2000,  # cell T104
            width_in=114.66,  # cell X104
            width_out=114.66,  # cell Y104
            height_in=111.27988185917508,  # cell Z104
            height_out=111.27988185917508,  # cell AA104
            mtop=2,  # cell AB104
            mbottom=2,  # cell AC104
            mleft=2,  # cell AD104
            mright=2,  # cell AE104
        ),
        "guide_06_28": guide(  # row 106
            x=82395.33178605184,  # cell H106
            y=702.3622679402808,  # cell I106
            z=-422.42085644514486,  # cell J106
            chi=0,  # cell F106
            phi=-0.7448451336700711,  # cell G106
            length=2000,  # cell T106
            width_in=114.66,  # cell X106
            width_out=114.66,  # cell Y106
            height_in=111.27988185917508,  # cell Z106
            height_out=111.27988185917508,  # cell AA106
            mtop=2,  # cell AB106
            mbottom=2,  # cell AC106
            mleft=2,  # cell AD106
            mright=2,  # cell AE106
        ),
        "guide_06_29": guide(  # row 108
            x=84397.1626194343,  # cell H108
            y=728.3875348808086,  # cell I108
            z=-422.42085644514486,  # cell J108
            chi=0,  # cell F108
            phi=-0.7448451336700711,  # cell G108
            length=2000,  # cell T108
            width_in=114.66,  # cell X108
            width_out=114.66,  # cell Y108
            height_in=111.27988185917508,  # cell Z108
            height_out=111.27988185917508,  # cell AA108
            mtop=2,  # cell AB108
            mbottom=2,  # cell AC108
            mleft=2,  # cell AD108
            mright=2,  # cell AE108
        ),
        "guide_06_30": guide(  # row 110
            x=86397.99353731556,  # cell H110
            y=754.3998021875,  # cell I110
            z=-422.42085644514486,  # cell J110
            chi=0,  # cell F110
            phi=-0.7448451336700711,  # cell G110
            length=2000,  # cell T110
            width_in=114.66,  # cell X110
            width_out=114.66,  # cell Y110
            height_in=111.27988185917508,  # cell Z110
            height_out=111.27988185917508,  # cell AA110
            mtop=2,  # cell AB110
            mbottom=2,  # cell AC110
            mleft=2,  # cell AD110
            mright=2,  # cell AE110
        ),
        "guide_06_31": guide(  # row 112
            x=88399.82437069801,  # cell H112
            y=780.4250691280278,  # cell I112
            z=-422.42085644514486,  # cell J112
            chi=0,  # cell F112
            phi=-0.7448451336700711,  # cell G112
            length=2000,  # cell T112
            width_in=114.66,  # cell X112
            width_out=114.66,  # cell Y112
            height_in=111.27988185917508,  # cell Z112
            height_out=111.27988185917508,  # cell AA112
            mtop=2,  # cell AB112
            mbottom=2,  # cell AC112
            mleft=2,  # cell AD112
            mright=2,  # cell AE112
        ),
        "guide_06_32": guide(  # row 114
            x=90400.65528857928,  # cell H114
            y=806.4373364347191,  # cell I114
            z=-422.42085644514486,  # cell J114
            chi=0,  # cell F114
            phi=-0.7448451336700711,  # cell G114
            length=2000,  # cell T114
            width_in=114.66,  # cell X114
            width_out=114.66,  # cell Y114
            height_in=111.27988185917508,  # cell Z114
            height_out=111.27988185917508,  # cell AA114
            mtop=2,  # cell AB114
            mbottom=2,  # cell AC114
            mleft=2,  # cell AD114
            mright=2,  # cell AE114
        ),
        "guide_06_331": guide(  # row 116
            x=92402.48612196173,  # cell H116
            y=832.4626033752469,  # cell I116
            z=-422.42085644514486,  # cell J116
            chi=0,  # cell F116
            phi=-0.7448451336700711,  # cell G116
            length=1162.36,  # cell T116
            width_in=114.66,  # cell X116
            width_out=114.66,  # cell Y116
            height_in=111.27988185917508,  # cell Z116
            height_out=111.27988185917508,  # cell AA116
            mtop=2,  # cell AB116
            mbottom=2,  # cell AC116
            mleft=2,  # cell AD116
            mright=2,  # cell AE116
        ),
        "guide_06_332": guide(  # row 117, converging in height, to row 129
            x=93564.74790392497,  # cell H117
            y=847.5728577613567,  # cell I117
            z=-422.42085644514486,  # cell J117
            chi=0,  # cell F117
            phi=-0.7448451336700711,  # cell G117
            length=837.6400000000001,  # cell T117
            width_in=114.66,  # cell X117
            width_out=114.66,  # cell Y117
            height_in=111.27569760764854,  # cell Z117
            height_out=109.39139494316318,  # cell AA117
            mtop=2.5,  # cell AB117
            mbottom=2.5,  # cell AC117
            mleft=2.5,  # cell AD117
            mright=2.5,  # cell AE117
        ),
        "guide_06_34": guide(  # row 119
            x=94403.31703984298,  # cell H119
            y=858.4748706819382,  # cell I119
            z=-422.42085644514486,  # cell J119
            chi=0,  # cell F119
            phi=-0.7448451336700711,  # cell G119
            length=2000,  # cell T119
            width_in=114.66,  # cell X119
            width_out=114.66,  # cell Y119
            height_in=109.3890760127605,  # cell Z119
            height_out=104.3987616196667,  # cell AA119
            mtop=2.5,  # cell AB119
            mbottom=2.5,  # cell AC119
            mleft=2.5,  # cell AD119
            mright=2.5,  # cell AE119
        ),
        "guide_06_35": guide(  # row 121
            x=96405.14787322543,  # cell H121
            y=884.500137622466,  # cell I121
            z=-422.42085644514486,  # cell J121
            chi=0,  # cell F121
            phi=-0.7448451336700711,  # cell G121
            length=2000,  # cell T121
            width_in=114.66,  # cell X121
            width_out=114.66,  # cell Y121
            height_in=104.39340179321815,  # cell Z121
            height_out=98.62388829152566,  # cell AA121
            mtop=2.5,  # cell AB121
            mbottom=2.5,  # cell AC121
            mleft=2.5,  # cell AD121
            mright=2.5,  # cell AE121
        ),
        "guide_06_36": guide(  # row 123
            x=98405.9787911067,  # cell H123
            y=910.5124049291574,  # cell I123
            z=-422.42085644514486,  # cell J123
            chi=0,  # cell F123
            phi=-0.7448451336700711,  # cell G123
            length=2000,  # cell T123
            width_in=114.66,  # cell X123
            width_out=114.66,  # cell Y123
            height_in=98.62078679804394,  # cell Z123
            height_out=91.92583036609614,  # cell AA123
            mtop=2.5,  # cell AB123
            mbottom=2.5,  # cell AC123
            mleft=2.5,  # cell AD123
            mright=2.5,  # cell AE123
        ),
        "guide_06_37": guide(  # row 125
            x=100407.80962448915,  # cell H125
            y=936.5376718696851,  # cell I125
            z=-422.42085644514486,  # cell J125
            chi=0,  # cell F125
            phi=-0.7448451336700711,  # cell G125
            length=1677.0349999999999,  # cell T125
            width_in=114.66,  # cell X125
            width_out=114.66,  # cell Y125
            height_in=91.91860728905549,  # cell Z125
            height_out=85.43312975229115,  # cell AA125
            mtop=2.5,  # cell AB125
            mbottom=2.5,  # cell AC125
            mleft=2.5,  # cell AD125
            mright=2.5,  # cell AE125
        ),
        "guide_06_38": guide(  # row 127
            x=102085.70283252856,  # cell H127
            y=958.3515124343946,  # cell I127
            z=-422.42085644514486,  # cell J127
            chi=0,  # cell F127
            phi=-0.7448451336700711,  # cell G127
            length=1677.0349999999999,  # cell T127
            width_in=114.66,  # cell X127
            width_out=114.66,  # cell Y127
            height_in=85.42898755223936,  # cell Z127
            height_out=77.93905381840912,  # cell AA127
            mtop=2.5,  # cell AB127
            mbottom=2.5,  # cell AC127
            mleft=2.5,  # cell AD127
            mright=2.5,  # cell AE127
        ),
        "guide_06_39": guide(  # row 129
            x=103764.59595606917,  # cell H129
            y=980.1783526329406,  # cell I129
            z=-422.42085644514486,  # cell J129
            chi=0,  # cell F129
            phi=-0.7448451336700711,  # cell G129
            length=800,  # cell T129
            width_in=114.66,  # cell X129
            width_out=114.66,  # cell Y129
            height_in=77.92941070199981,  # cell Z129
            height_out=73.91482959988588,  # cell AA129
            mtop=2.5,  # cell AB129
            mbottom=2.5,  # cell AC129
            mleft=2.5,  # cell AD129
            mright=2.5,  # cell AE129
        ),
        "arm_bwc3": arm(  # row 130, the gap of bandwidth chopper 3
            x=104564.53,  # cell H130
            y=990.5780597020826,  # cell I130
            z=-422.42085644514486,  # cell J130
            chi=0,  # cell F130
            phi=-0.7448451336700711,  # cell G130
        ),
    }
    to_psc: dict[str, Component] = {
        "guide_06_40": guide(  # row 131, constant, to row 137
            x=104594.5274650357,  # cell H131
            y=990.9680487171754,  # cell I131
            z=-422.42085644514486,  # cell J131
            chi=0,  # cell F131
            phi=-0.7448451336700711,  # cell G131
            length=500,  # cell T131
            width_in=114.66,  # cell X131
            width_out=114.66,  # cell Y131
            height_in=73.89,  # cell Z131
            height_out=73.89,  # cell AA131
            mtop=2.5,  # cell AB131
            mbottom=2.5,  # cell AC131
            mleft=2.5,  # cell AD131
            mright=2.5,  # cell AE131
        ),
        "guide_07_01": guide(  # row 133
            x=105104.48437064262,  # cell H133
            y=997.5978619737534,  # cell I133
            z=-422.42085644514486,  # cell J133
            chi=0,  # cell F133
            phi=-0.7448451336700711,  # cell G133
            length=518.88,  # cell T133
            width_in=114.66,  # cell X133
            width_out=114.66,  # cell Y133
            height_in=73.89,  # cell Z133
            height_out=73.89,  # cell AA133
            mtop=2.5,  # cell AB133
            mbottom=2.5,  # cell AC133
            mleft=2.5,  # cell AD133
            mright=2.5,  # cell AE133
        ),
        "arm_psc": arm(  # row 134, the gap of the pulse-shaping choppers
            x=105623.31999999999,  # cell H134
            y=1004.3431119787989,  # cell I134
            z=-422.42085644514486,  # cell J134
            chi=0,  # cell F134
            phi=-0.7448451336700711,  # cell G134
        ),
    }
    to_mc: dict[str, Component] = {
        "guide_08_01": guide(  # row 135
            x=105730.62093243269,  # cell H135
            y=1005.738102685786,  # cell I135
            z=-422.42085644514486,  # cell J135
            chi=0,  # cell F135
            phi=-0.7448451336700711,  # cell G135
            length=2480.2,  # cell T135
            width_in=114.66,  # cell X135
            width_out=114.66,  # cell Y135
            height_in=73.89,  # cell Z135
            height_out=73.89,  # cell AA135
            mtop=2.5,  # cell AB135
            mbottom=2.5,  # cell AC135
            mleft=2,  # cell AD135
            mright=2,  # cell AE135
        ),
        "guide_08_02": guide(  # row 137
            x=108211.6112739854,  # cell H137
            y=1037.9927941607298,  # cell I137
            z=-422.42085644514486,  # cell J137
            chi=0,  # cell F137
            phi=-0.7448451336700711,  # cell G137
            length=2480.2,  # cell T137
            width_in=114.66,  # cell X137
            width_out=114.46,  # cell Y137
            height_in=73.89,  # cell Z137
            height_out=73.89,  # cell AA137
            mtop=2.5,  # cell AB137
            mbottom=2.5,  # cell AC137
            mleft=2,  # cell AD137
            mright=2,  # cell AE137
        ),
        "guide_08_03": guide(  # row 139, converging in width, and in height from row 153
            x=110693.6015310393,  # cell H139
            y=1070.26048526951,  # cell I139
            z=-422.42085644514486,  # cell J139
            chi=0,  # cell F139
            phi=-0.7448451336700711,  # cell G139
            length=2000,  # cell T139
            width_in=114.46,  # cell X139
            width_out=114.36332854461857,  # cell Y139
            height_in=73.89,  # cell Z139
            height_out=73.89,  # cell AA139
            mtop=3,  # cell AB139
            mbottom=3,  # cell AC139
            mleft=2,  # cell AD139
            mright=2,  # cell AE139
        ),
        "guide_08_04": guide(  # row 141
            x=112694.43244892056,  # cell H141
            y=1096.2727525762011,  # cell I141
            z=-422.42085644514486,  # cell J141
            chi=0,  # cell F141
            phi=-0.7448451336700711,  # cell G141
            length=2000,  # cell T141
            width_in=114.36323180808596,  # cell X141
            width_out=114.07262869031695,  # cell Y141
            height_in=73.89,  # cell Z141
            height_out=73.89,  # cell AA141
            mtop=3,  # cell AB141
            mbottom=3,  # cell AC141
            mleft=2,  # cell AD141
            mright=2,  # cell AE141
        ),
        "guide_08_05": guide(  # row 143
            x=114696.26328230301,  # cell H143
            y=1122.2980195167288,  # cell I143
            z=-422.42085644514486,  # cell J143
            chi=0,  # cell F143
            phi=-0.7448451336700711,  # cell G143
            length=2000,  # cell T143
            width_in=114.07224066064788,  # cell X143
            width_out=113.5861184578029,  # cell Y143
            height_in=73.89,  # cell Z143
            height_out=73.89,  # cell AA143
            mtop=3,  # cell AB143
            mbottom=3,  # cell AC143
            mleft=2,  # cell AD143
            mright=2,  # cell AE143
        ),
        "guide_08_06": guide(  # row 145
            x=116697.09420018428,  # cell H145
            y=1148.31028682342,  # cell I145
            z=-422.42085644514486,  # cell J145
            chi=0,  # cell F145
            phi=-0.7448451336700711,  # cell G145
            length=2000,  # cell T145
            width_in=113.58582616482678,  # cell X145
            width_out=112.90175655604327,  # cell Y145
            height_in=73.89,  # cell Z145
            height_out=73.89,  # cell AA145
            mtop=3,  # cell AB145
            mbottom=3,  # cell AC145
            mleft=2,  # cell AD145
            mright=2,  # cell AE145
        ),
        "guide_08_07": guide(  # row 147
            x=118698.92503356673,  # cell H147
            y=1174.3355537639477,  # cell I147
            z=-422.42085644514486,  # cell J147
            chi=0,  # cell F147
            phi=-0.7448451336700711,  # cell G147
            length=2000,  # cell T147
            width_in=112.90097234905937,  # cell X147
            width_out=112.01522519504762,  # cell Y147
            height_in=73.89,  # cell Z147
            height_out=73.89,  # cell AA147
            mtop=3,  # cell AB147
            mbottom=3,  # cell AC147
            mleft=2,  # cell AD147
            mright=2,  # cell AE147
        ),
        "guide_08_08": guide(  # row 149
            x=120699.755951448,  # cell H149
            y=1200.347821070639,  # cell I149
            z=-422.42085644514486,  # cell J149
            chi=0,  # cell F149
            phi=-0.7448451336700711,  # cell G149
            length=2000,  # cell T149
            width_in=112.01473117489516,  # cell X149
            width_out=110.92257493063678,  # cell Y149
            height_in=73.89,  # cell Z149
            height_out=73.89,  # cell AA149
            mtop=3,  # cell AB149
            mbottom=3,  # cell AC149
            mleft=2,  # cell AD149
            mright=2,  # cell AE149
        ),
        "guide_08_09": guide(  # row 151
            x=122701.58678483045,  # cell H151
            y=1226.3730880111666,  # cell I151
            z=-422.42085644514486,  # cell J151
            chi=0,  # cell F151
            phi=-0.7448451336700711,  # cell G151
            length=2000,  # cell T151
            width_in=110.92137757919757,  # cell X151
            width_out=109.61653170651755,  # cell Y151
            height_in=73.89,  # cell Z151
            height_out=73.89,  # cell AA151
            mtop=3,  # cell AB151
            mbottom=3,  # cell AC151
            mleft=2,  # cell AD151
            mright=2,  # cell AE151
        ),
        "guide_08_10": guide(  # row 153
            x=124702.41770271171,  # cell H153
            y=1252.3853553178578,  # cell I153
            z=-422.42085644514486,  # cell J153
            chi=0,  # cell F153
            phi=-0.7448451336700711,  # cell G153
            length=2000,  # cell T153
            width_in=109.61582492278816,  # cell X153
            width_out=108.09069156836284,  # cell Y153
            height_in=73.88999999888553,  # cell Z153
            height_out=73.81710248034494,  # cell AA153
            mtop=3,  # cell AB153
            mbottom=3,  # cell AC153
            mleft=2,  # cell AD153
            mright=2,  # cell AE153
        ),
        "guide_08_11": guide(  # row 155
            x=126704.24853609416,  # cell H155
            y=1278.4106222583855,  # cell I155
            z=-422.42085644514486,  # cell J155
            chi=0,  # cell F155
            phi=-0.7448451336700711,  # cell G155
            length=2000,  # cell T155
            width_in=108.08905323598819,  # cell X155
            width_out=106.33403194610952,  # cell Y155
            height_in=73.81695652218018,  # cell Z155
            height_out=73.59764778047658,  # cell AA155
            mtop=3,  # cell AB155
            mbottom=3,  # cell AC155
            mleft=2,  # cell AD155
            mright=2,  # cell AE155
        ),
        "guide_08_12": guide(  # row 157
            x=128705.07945397543,  # cell H157
            y=1304.4228895650767,  # cell I157
            z=-422.42085644514486,  # cell J157
            chi=0,  # cell F157
            phi=-0.7448451336700711,  # cell G157
            length=2000,  # cell T157
            width_in=106.3330951562316,  # cell X157
            width_out=104.33669803522007,  # cell Y157
            height_in=73.59750135975,  # cell Z157
            height_out=73.23053895082423,  # cell AA157
            mtop=3.5,  # cell AB157
            mbottom=3.5,  # cell AC157
            mleft=2,  # cell AD157
            mright=2,  # cell AE157
        ),
        "guide_08_13": guide(  # row 159
            x=130706.91028735788,  # cell H159
            y=1330.4481565056044,  # cell I159
            z=-422.42085644514486,  # cell J159
            chi=0,  # cell F159
            phi=-0.7448451336700711,  # cell G159
            length=2000,  # cell T159
            width_in=104.33457640456113,  # cell X159
            width_out=102.08250415790185,  # cell Y159
            height_in=73.2300974582292,  # cell Z159
            height_out=72.71316918024365,  # cell AA159
            mtop=3.5,  # cell AB159
            mbottom=3.5,  # cell AC159
            mleft=2,  # cell AD159
            mright=2,  # cell AE159
        ),
        "guide_08_14": guide(  # row 161
            x=132707.74120523914,  # cell H161
            y=1356.4604238122956,  # cell I161
            z=-422.42085644514486,  # cell J161
            chi=0,  # cell F161
            phi=-0.7448451336700711,  # cell G161
            length=2000,  # cell T161
            width_in=102.08131149327282,  # cell X161
            width_out=99.55633826502452,  # cell Y161
            height_in=72.71287274870237,  # cell Z161
            height_out=72.04282477714283,  # cell AA161
            mtop=3.5,  # cell AB161
            mbottom=3.5,  # cell AC161
            mleft=2,  # cell AD161
            mright=2,  # cell AE161
        ),
        "guide_08_15": guide(  # row 163
            x=134709.57203862158,  # cell H163
            y=1382.4856907528233,  # cell I163
            z=-422.42085644514486,  # cell J163
            chi=0,  # cell F163
            phi=-0.7448451336700711,  # cell G163
            length=2000,  # cell T163
            width_in=99.55367002912378,  # cell X163
            width_out=96.73426455533428,  # cell Y163
            height_in=72.0420767863721,  # cell Z163
            height_out=71.21450509592736,  # cell AA163
            mtop=3.5,  # cell AB163
            mbottom=3.5,  # cell AC163
            mleft=2,  # cell AD163
            mright=2,  # cell AE163
        ),
        "guide_08_16": guide(  # row 165
            x=136710.40295650283,  # cell H165
            y=1408.4979580595145,  # cell I165
            z=-422.42085644514486,  # cell J165
            chi=0,  # cell F165
            phi=-0.7448451336700711,  # cell G165
            length=2000,  # cell T165
            width_in=96.73277710037343,  # cell X165
            width_out=93.59247512096694,  # cell Y165
            height_in=71.21405107664532,  # cell Z165
            height_out=70.22346416347962,  # cell AA165
            mtop=3.5,  # cell AB165
            mbottom=3.5,  # cell AC165
            mleft=2.5,  # cell AD165
            mright=2.5,  # cell AE165
        ),
        "guide_08_17": guide(  # row 167
            x=138712.23378988527,  # cell H167
            y=1434.5232250000422,  # cell I167
            z=-422.42085644514486,  # cell J167
            chi=0,  # cell F167
            phi=-0.7448451336700711,  # cell G167
            length=2000,  # cell T167
            width_in=93.58916378351532,  # cell X167
            width_out=90.0942129199001,  # cell Y167
            height_in=70.2223898201479,  # cell Z167
            height_out=69.06168293843176,  # cell AA167
            mtop=3.5,  # cell AB167
            mbottom=3.5,  # cell AC167
            mleft=2.5,  # cell AD167
            mright=2.5,  # cell AE167
        ),
        "guide_08_18": guide(  # row 169
            x=140713.06470776652,  # cell H169
            y=1460.5354923067334,  # cell I169
            z=-422.42085644514486,  # cell J169
            chi=0,  # cell F169
            phi=-0.7448451336700711,  # cell G169
            length=2000,  # cell T169
            width_in=90.09237011865103,  # cell X169
            width_out=86.19981025779111,  # cell Y169
            height_in=69.0610586983185,  # cell Z169
            height_out=67.72156810630119,  # cell AA169
            mtop=3.5,  # cell AB169
            mbottom=3.5,  # cell AC169
            mleft=2.5,  # cell AD169
            mright=2.5,  # cell AE169
        ),
        "guide_08_19": guide(  # row 171
            x=142714.89554114896,  # cell H171
            y=1486.560759247261,  # cell I171
            z=-422.42085644514486,  # cell J171
            chi=0,  # cell F171
            phi=-0.7448451336700711,  # cell G171
            length=2000,  # cell T171
            width_in=86.19570127221654,  # cell X171
            width_out=81.84855186826901,  # cell Y171
            height_in=67.72013575606273,  # cell Z171
            height_out=66.19090507710686,  # cell AA171
            mtop=3.5,  # cell AB171
            mbottom=3.5,  # cell AC171
            mleft=3,  # cell AD171
            mright=3,  # cell AE171
        ),
        "guide_08_20": guide(  # row 173
            x=144715.7264590302,  # cell H173
            y=1512.5730265539523,  # cell I173
            z=-422.42085644514486,  # cell J173
            chi=0,  # cell F173
            phi=-0.7448451336700711,  # cell G173
            length=2000,  # cell T173
            width_in=81.84625294002784,  # cell X173
            width_out=76.96773451499725,  # cell Y173
            height_in=66.1900909254975,  # cell Z173
            height_out=64.45770865335473,  # cell AA173
            mtop=3.5,  # cell AB173
            mbottom=3.5,  # cell AC173
            mleft=3,  # cell AD173
            mright=3,  # cell AE173
        ),
        "guide_08_21": guide(  # row 175
            x=146717.55729241265,  # cell H175
            y=1538.59829349448,  # cell I175
            z=-422.42085644514486,  # cell J175
            chi=0,  # cell F175
            phi=-0.7448451336700711,  # cell G175
            length=2000,  # cell T175
            width_in=76.9625573678673,  # cell X175
            width_out=71.44349279599015,  # cell Y175
            height_in=64.4558693370418,  # cell Z175
            height_out=62.503322551462404,  # cell AA175
            mtop=3.5,  # cell AB175
            mbottom=3.5,  # cell AC175
            mleft=3,  # cell AD175
            mright=3,  # cell AE175
        ),
        "guide_08_22": guide(  # row 177
            x=148718.3882102939,  # cell H177
            y=1564.6105608011712,  # cell I177
            z=-422.42085644514486,  # cell J177
            chi=0,  # cell F177
            phi=-0.7448451336700711,  # cell G177
            length=2000,  # cell T177
            width_in=71.44054917014267,  # cell X177
            width_out=65.11855344634543,  # cell Y177
            height_in=62.50228792073246,  # cell Z177
            height_out=60.30830103155461,  # cell AA177
            mtop=3.5,  # cell AB177
            mbottom=3.5,  # cell AC177
            mleft=3.5,  # cell AD177
            mright=3.5,  # cell AE177
        ),
        "guide_08_23": guide(  # row 179
            x=150720.21904367633,  # cell H179
            y=1590.635827741699,  # cell I179
            z=-422.42085644514486,  # cell J179
            chi=0,  # cell F179
            phi=-0.7448451336700711,  # cell G179
            length=2000,  # cell T179
            width_in=65.11175419290193,  # cell X179
            width_out=57.72286255488553,  # cell Y179
            height_in=60.30597771013691,  # cell Z179
            height_out=57.84292779091675,  # cell AA179
            mtop=3.5,  # cell AB179
            mbottom=3.5,  # cell AC179
            mleft=3.5,  # cell AD179
            mright=3.5,  # cell AE179
        ),
        "guide_08_24": guide(  # row 181
            x=152721.04996155758,  # cell H181
            y=1616.6480950483901,  # cell I181
            z=-422.42085644514486,  # cell J181
            chi=0,  # cell F181
            phi=-0.7448451336700711,  # cell G181
            length=2000,  # cell T181
            width_in=57.71883565892217,  # cell X181
            width_out=48.78102938847784,  # cell Y181
            height_in=57.841623458666305,  # cell Z181
            height_out=55.07354802617992,  # cell AA181
            mtop=3.5,  # cell AB181
            mbottom=3.5,  # cell AC181
            mleft=3.5,  # cell AD181
            mright=3.5,  # cell AE181
        ),
        "guide_08_25": guide(  # row 183
            x=154722.88079494002,  # cell H183
            y=1642.6733619889178,  # cell I183
            z=-422.42085644514486,  # cell J183
            chi=0,  # cell F183
            phi=-0.7448451336700711,  # cell G183
            length=1604.4,  # cell T183
            width_in=48.77104479333032,  # cell X183
            width_out=39.78699087056744,  # cell Y183
            height_in=55.070612431357034,  # cell Z183
            height_out=52.5970589597697,  # cell AA183
            mtop=3.5,  # cell AB183
            mbottom=3.5,  # cell AC183
            mleft=3.5,  # cell AD183
            mright=3.5,  # cell AE183
        ),
        "guide_08_26": guide(  # row 185
            x=156328.1451405505,  # cell H185
            y=1663.5429741499183,  # cell I185
            z=-422.42085644514486,  # cell J185
            chi=0,  # cell F185
            phi=-0.7448451336700711,  # cell G185
            length=1604.4,  # cell T185
            width_in=39.78064694984217,  # cell X185
            width_out=27.539091112840772,  # cell Y185
            height_in=52.5954398595165,  # cell Z185
            height_out=49.8606879242676,  # cell AA185
            mtop=3.5,  # cell AB185
            mbottom=3.5,  # cell AC185
            mleft=3.5,  # cell AD185
            mright=3.5,  # cell AE185
        ),
        "guide_09_01": guide(  # row 187
            x=157942.4087256717,  # cell H187
            y=1684.5295830154466,  # cell I187
            z=-422.42085644514486,  # cell J187
            chi=0,  # cell F187
            phi=-0.7448451336700711,  # cell G187
            length=476.3,  # cell T187
            width_in=27.4440480069477,  # cell X187
            width_out=22.42747203925759,  # cell Y187
            height_in=49.84273602025601,  # cell Z187
            height_out=48.97370053838532,  # cell AA187
            mtop=3.5,  # cell AB187
            mbottom=3.5,  # cell AC187
            mleft=3.5,  # cell AD187
            mright=3.5,  # cell AE187
        ),
        "arm_mc": arm(  # row 188, the gap of the monochromatic choppers
            x=158418.66845151584,  # cell H188
            y=1690.721308611737,  # cell I188
            z=-422.42085644514486,  # cell J188
            chi=0,  # cell F188
            phi=-0.7448451336700711,  # cell G188
        ),
    }
    to_guide_exchanger: dict[str, Component] = {
        "guide_09_02": guide(  # row 189
            x=158491.63228563766,  # cell H189
            y=1691.6698918927812,  # cell I189
            z=-422.42085644514486,  # cell J189
            chi=0,  # cell F189
            phi=-0.7448451336700711,  # cell G189
            length=300,  # cell T189
            width_in=22.58,  # cell X189
            width_out=22.58,  # cell Y189
            height_in=49,  # cell Z189
            height_out=49,  # cell AA189
            mtop=3.5,  # cell AB189
            mbottom=3.5,  # cell AC189
            mleft=3.5,  # cell AD189
            mright=3.5,  # cell AE189
        ),
    }
    return [to_bwc3, to_psc, to_mc, to_guide_exchanger]


def sample_pos() -> Arm:
    """Create the arm at the sample position.

    Returns
    -------
    :
        A new arm, relative to ``arm_iscs``, along the trajectory of the last guide.
    """
    return arm(x=SAMPLE_X, y=SAMPLE_Y, z=SAMPLE_Z, phi=GUIDE_09_02_PHI)
