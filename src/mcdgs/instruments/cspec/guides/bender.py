"""The guides of CSPEC from the NBOA to the end of the S-bender, rows 14 to 76.

Dimensions are from the Specification sheet of CSPEC_nl-specs_1.E.xlsx, in ``reference/``, whose
rows and cells are named in the comments, in mm and deg as the sheet gives them. The values are the
sheet's, not rounded. Each guide and gap is placed by its own start in ISCS: see :mod:`.frame`. The
gaps of 1 to 12 mm between the elements are left empty.

Over these rows the trajectory turns down and back up in a vertical S-bender, 1600 m in radius,
which takes the beam 422 mm down, and sideways in a horizontal curve, 4000 m in radius, by
0.745 deg. Both end at row 77, after which the guide is straight: see :mod:`.straight`.

The gaps of the bandwidth choppers 1 and 2 are arms, for the choppers to be placed relative to.
"""

from mcdgs.core.components import Component

from .frame import arm, guide

# the start of the guide system in the McStas frame of the source, the entrance of W03-01-011:
# the source focuses on it
GUIDE_START_X = 0.0163  # m, cell B7 of Definitions
GUIDE_START_Z = 1.90424  # m, cell D7 of Definitions
GUIDE_01_011_IN_WIDTH = 45  # mm, cell X14
GUIDE_01_011_IN_HEIGHT = 43  # mm, cell Z14

# W03-01-01 is one element in two rows, 14 and 15, and the sheet gives only its entry and its
# exit. Its middle is where a linear taper puts it. FIXME, not in CSPEC_nl-specs_1.E.xlsx
GUIDE_01_012_OUT_WIDTH = 61.7  # mm, cell Y15
GUIDE_01_012_OUT_HEIGHT = 51  # mm, cell AA15
GUIDE_01_01_MIDDLE_WIDTH = (GUIDE_01_011_IN_WIDTH + GUIDE_01_012_OUT_WIDTH) / 2
GUIDE_01_01_MIDDLE_HEIGHT = (GUIDE_01_011_IN_HEIGHT + GUIDE_01_012_OUT_HEIGHT) / 2

# the flight path from the moderator to the start of each chopper gap, along the trajectory
ARM_BWC1_PATH = 14.918415010646674  # m, cell U31
ARM_BWC2_PATH = 20.440446374181604  # m, cell U39


def guides_bender() -> list[dict[str, Component]]:
    """Create the guides from the NBOA to the end of the S-bender.

    Returns
    -------
    :
        Components by name, in beam order, in three pieces: up to ``arm_bwc1``, up to
        ``arm_bwc2``, and the rest.
    """
    to_bwc1: dict[str, Component] = {
        "guide_01_011": guide(  # row 14, straight: the trajectory turns from row 15 on
            x=1879.938878035614,  # cell H14
            y=-0.0050296957369937445,  # cell I14
            z=0,  # cell J14
            chi=0,  # cell F14
            phi=0,  # cell G14
            length=499.5,  # cell T14
            width_in=GUIDE_01_011_IN_WIDTH,
            width_out=GUIDE_01_01_MIDDLE_WIDTH,
            height_in=GUIDE_01_011_IN_HEIGHT,
            height_out=GUIDE_01_01_MIDDLE_HEIGHT,
            mtop=4,  # cell AB14
            mbottom=3.5,  # cell AC14
            mleft=3,  # cell AD14
            mright=3.5,  # cell AE14
        ),
        "guide_01_012": guide(  # row 15
            x=2379.438878035614,  # cell H15
            y=-0.0050296957369937445,  # cell I15
            z=0,  # cell J15
            chi=-0.008943513083370195,  # cell F15
            phi=-0.0035774052333480775,  # cell G15
            length=499.5,  # cell T15
            width_in=GUIDE_01_01_MIDDLE_WIDTH,
            width_out=GUIDE_01_012_OUT_WIDTH,
            height_in=GUIDE_01_01_MIDDLE_HEIGHT,
            height_out=GUIDE_01_012_OUT_HEIGHT,
            mtop=4,  # cell AB15
            mbottom=3.5,  # cell AC15
            mleft=3,  # cell AD15
            mright=3.5,  # cell AE15
        ),
        "guide_01_02": guide(  # row 17
            x=2879.9388709201144,  # cell H17
            y=0.02628283549241694,  # cell I17
            z=-0.07828132780329196,  # cell J17
            chi=-0.03580986219567646,  # cell F17
            phi=-0.014323944878270581,  # cell G17
            length=999,  # cell T17
            width_in=61.7,  # cell X17
            width_out=70,  # cell Y17
            height_in=51,  # cell Z17
            height_out=55,  # cell AA17
            mtop=4,  # cell AB17
            mbottom=3.5,  # cell AC17
            mleft=3,  # cell AD17
            mright=3.5,  # cell AE17
        ),
        "guide_01_03": guide(  # row 19
            x=3879.9386440744242,  # cell H19
            y=0.2764078328820654,  # cell I19
            z=-0.7035937870165497,  # cell J19
            chi=-0.0802319962494131,  # cell F19
            phi=-0.032092798499765235,  # cell G19
            length=1480,  # cell T19
            width_in=70,  # cell X19
            width_out=70,  # cell Y19
            height_in=55,  # cell Z19
            height_out=55,  # cell AA19
            mtop=3.5,  # cell AB19
            mbottom=3,  # cell AC19
            mleft=3,  # cell AD19
            mright=3.5,  # cell AE19
        ),
        "guide_02_01": guide(  # row 21, BBG
            x=5361.93695683173,  # cell H21
            y=1.1068835393963987,  # cell I21
            z=-2.779782482553066,  # cell J21
            chi=-0.11604185844508956,  # cell F21
            phi=-0.04641674337803581,  # cell G21
            length=516,  # cell T21
            width_in=70,  # cell X21
            width_out=70,  # cell Y21
            height_in=55,  # cell Z21
            height_out=55,  # cell AA21
            mtop=3.5,  # cell AB21
            mbottom=3,  # cell AC21
            mleft=3,  # cell AD21
            mright=3.5,  # cell AE21
        ),
        "guide_03_01": guide(  # row 23
            x=5902.935659397056,  # cell H23
            y=1.5468517408535696,  # cell I23
            z=-3.8797023490714917,  # cell J23
            chi=-0.16093334493477746,  # cell F23
            phi=-0.06437333797391095,  # cell G23
            length=1941.2135851507358,  # cell T23
            width_in=70,  # cell X23
            width_out=70,  # cell Y23
            height_in=55,  # cell Z23
            height_out=55,  # cell AA23
            mtop=3.5,  # cell AB23
            mbottom=3,  # cell AC23
            mleft=3,  # cell AD23
            mright=3.5,  # cell AE23
        ),
        "guide_03_02": guide(  # row 25
            x=7845.140355026644,  # cell H25
            y=3.72922287468955,  # cell I25
            z=-9.335624155657364,  # cell J25
            chi=-0.23048374577359604,  # cell F25
            phi=-0.09219349830943838,  # cell G25
            length=1941.2135851507358,  # cell T25
            width_in=70,  # cell X25
            width_out=70,  # cell Y25
            height_in=55,  # cell Z25
            height_out=55,  # cell AA25
            mtop=3.5,  # cell AB25
            mbottom=3,  # cell AC25
            mleft=3,  # cell AD25
            mright=3.5,  # cell AE25
        ),
        "guide_04_01": guide(  # row 27
            x=9798.33557131201,  # cell H27
            y=6.875028434485122,  # cell I27
            z=-17.20012019700441,  # cell J27
            chi=-0.30042805509656706,  # cell F27
            phi=-0.12017122203862676,  # cell G27
            length=1941.2135851507358,  # cell T27
            width_in=70,  # cell X27
            width_out=70,  # cell Y27
            height_in=55,  # cell Z27
            height_out=55,  # cell AA27
            mtop=3.5,  # cell AB27
            mbottom=3,  # cell AC27
            mleft=3,  # cell AD27
            mright=3.5,  # cell AE27
        ),
        "guide_04_02": guide(  # row 29
            x=11740.518181289499,  # cell H29
            y=10.948834932916247,  # cell I29
            z=-27.384597235968933,  # cell J29
            chi=-0.36997845593538564,  # cell F29
            phi=-0.14799138237415418,  # cell G29
            length=1941.2135851507358,  # cell T29
            width_in=70,  # cell X29
            width_out=70,  # cell Y29
            height_in=55,  # cell Z29
            height_out=55,  # cell AA29
            mtop=3.5,  # cell AB29
            mbottom=3,  # cell AC29
            mleft=3,  # cell AD29
            mright=3.5,  # cell AE29
        ),
        "guide_05_01": guide(  # row 31
            x=13693.684471934866,  # cell H31
            y=15.996786550155495,  # cell I31
            z=-40.00440251111285,  # cell J31
            chi=-0.4266513870874511,  # cell F31
            phi=-0.17066055483498035,  # cell G31
            length=1200,  # cell T31
            width_in=70,  # cell X31
            width_out=70,  # cell Y31
            height_in=55,  # cell Z31
            height_out=55,  # cell AA31
            mtop=3.5,  # cell AB31
            mbottom=3,  # cell AC31
            mleft=3,  # cell AD31
            mright=3.5,  # cell AE31
        ),
        "arm_bwc1": arm(  # row 32, the gap of bandwidth chopper 1
            x=14893.65,  # cell H32
            y=19.571087567137116,  # cell I32
            z=-48.94008568551895,  # cell J32
            chi=-0.44867445233779213,  # cell F32
            phi=-0.17946978093511678,  # cell G32
        ),
    }
    to_bwc2: dict[str, Component] = {
        "guide_05_02": guide(  # row 33
            x=14923.648933003813,  # cell H33
            y=19.66505757102631,  # cell I33
            z=-49.17500867839066,  # cell J33
            chi=-0.46353554514899786,  # cell F33
            phi=-0.1854142180595991,  # cell G33
            length=800,  # cell T33
            width_in=70,  # cell X33
            width_out=70,  # cell Y33
            height_in=55,  # cell Z33
            height_out=55,  # cell AA33
            mtop=3.5,  # cell AB33
            mbottom=3,  # cell AC33
            mleft=3,  # cell AD33
            mright=3.5,  # cell AE33
        ),
        "guide_05_03": guide(  # row 35
            x=15725.618483006216,  # cell H35
            y=22.260596585381638,  # cell I35
            z=-55.66379674594441,  # cell J35
            chi=-0.5127296214665215,  # cell F35
            phi=-0.20509184858660853,  # cell G35
            length=1943.5155334980932,  # cell T35
            width_in=70,  # cell X35
            width_out=70,  # cell Y35
            height_in=55,  # cell Z35
            height_out=55,  # cell AA35
            mtop=3.5,  # cell AB35
            mbottom=3,  # cell AC35
            mleft=3,  # cell AD35
            mright=3.5,  # cell AE35
        ),
        "guide_05_04": guide(  # row 37
            x=17670.043693710762,  # cell H37
            y=29.221272634885164,  # cell I37
            z=-73.06529175795642,  # cell J37
            chi=-0.5823624547584404,  # cell F37
            phi=-0.23294498190337612,  # cell G37
            length=1943.5155334980932,  # cell T37
            width_in=70,  # cell X37
            width_out=70,  # cell Y37
            height_in=55,  # cell Z37
            height_out=55,  # cell AA37
            mtop=3.5,  # cell AB37
            mbottom=3,  # cell AC37
            mleft=3,  # cell AD37
            mright=3.5,  # cell AE37
        ),
        "guide_05_05": guide(  # row 39
            x=19615.442639499255,  # cell H39
            y=37.13153589129047,  # cell I39
            z=-92.84066383883619,  # cell J39
            chi=-0.631556531075964,  # cell F39
            phi=-0.25262261243038553,  # cell G39
            length=800,  # cell T39
            width_in=70,  # cell X39
            width_out=70,  # cell Y39
            height_in=55,  # cell Z39
            height_out=55,  # cell AA39
            mtop=3.5,  # cell AB39
            mbottom=3,  # cell AC39
            mleft=3,  # cell AD39
            mright=3.5,  # cell AE39
        ),
        "arm_bwc2": arm(  # row 40, the gap of bandwidth chopper 2
            x=20415.39,  # cell H40
            y=40.65880154436686,  # cell I40
            z=-101.6586779741183,  # cell J40
            chi=-0.6464176238871697,  # cell F40
            phi=-0.2585670495548678,  # cell G40
        ),
    }
    to_end_of_bender: dict[str, Component] = {
        "guide_05_06": guide(  # row 41
            x=20445.3877852614,  # cell H41
            y=40.79418647538609,  # cell I41
            z=-101.99713427027166,  # cell J41
            chi=-0.6612787166983755,  # cell F41
            phi=-0.2645114866793501,  # cell G41
            length=800,  # cell T41
            width_in=70,  # cell X41
            width_out=70,  # cell Y41
            height_in=55,  # cell Z41
            height_out=55,  # cell AA41
            mtop=3.5,  # cell AB41
            mbottom=3,  # cell AC41
            mleft=3,  # cell AD41
            mright=3.5,  # cell AE41
        ),
        "guide_05_07": guide(  # row 43
            x=21247.325817546916,  # cell H43
            y=44.496884095530326,  # cell I43
            z=-111.25370567449023,  # cell J43
            chi=-0.7016927798444726,  # cell F43
            phi=-0.280677111937789,  # cell G43
            length=1453.147091226541,  # cell T43
            width_in=70,  # cell X43
            width_out=70,  # cell Y43
            height_in=55,  # cell Z43
            height_out=55,  # cell AA43
            mtop=3.5,  # cell AB43
            mbottom=3,  # cell AC43
            mleft=3,  # cell AD43
            mright=3.5,  # cell AE43
        ),
        "guide_05_08": guide(  # row 45
            x=22701.34640648984,  # cell H45
            y=51.62052551165356,  # cell I45
            z=-129.06243524182656,  # cell J45
            chi=-0.7537655867935388,  # cell F45
            phi=-0.30150623471741556,  # cell G45
            length=1453.147091226541,  # cell T45
            width_in=70,  # cell X45
            width_out=70,  # cell Y45
            height_in=55,  # cell Z45
            height_out=55,  # cell AA45
            mtop=3.5,  # cell AB45
            mbottom=3,  # cell AC45
            mleft=3,  # cell AD45
            mright=3.5,  # cell AE45
        ),
        "guide_05_09": guide(  # row 47, BWFA, the bunker wall feedthrough
            x=24157.35,  # cell H47
            y=59.38,  # cell I47
            z=-148.15,  # cell J47
            chi=-0.8164175743631511,  # cell F47
            phi=-0.32656702974526053,  # cell G47
            length=2040,  # cell T47
            width_in=70,  # cell X47
            width_out=70,  # cell Y47
            height_in=55,  # cell Z47
            height_out=55,  # cell AA47
            mtop=3.5,  # cell AB47
            mbottom=3,  # cell AC47
            mleft=3,  # cell AD47
            mright=3.5,  # cell AE47
        ),
        "guide_05_10": guide(  # row 49
            x=26198.10964269857,  # cell H49
            y=71.01321839605293,  # cell I49
            z=-177.2322192577787,  # cell J49
            chi=-0.8825225799763698,  # cell F49
            phi=-0.35300903199054806,  # cell G49
            length=1650,  # cell T49
            width_in=70,  # cell X49
            width_out=70,  # cell Y49
            height_in=55,  # cell Z49
            height_out=55,  # cell AA49
            mtop=3.5,  # cell AB49
            mbottom=3,  # cell AC49
            mleft=3,  # cell AD49
            mright=3.5,  # cell AE49
        ),
        "guide_06_011": guide(  # row 51
            x=27887.876714137048,  # cell H51
            y=81.43397948136817,  # cell I51
            z=-203.28325518691872,  # cell J51
            chi=-0.922277263931609,  # cell F51
            phi=-0.36891090557264367,  # cell G51
            length=490.32040994779163,  # cell T51
            width_in=70,  # cell X51
            width_out=70,  # cell Y51
            height_in=55,  # cell Z51
            height_out=55,  # cell AA51
            mtop=3.5,  # cell AB51
            mbottom=3,  # cell AC51
            mleft=3,  # cell AD51
            mright=3.5,  # cell AE51
        ),
        "guide_06_012": guide(  # row 52, the S-bender turns back up from here
            x=28378.12344068548,  # cell H52
            y=84.59098856952507,  # cell I52
            z=-211.17549160945825,  # cell J52
            chi=-0.8955688436516727,  # cell F52
            phi=-0.3866175962094014,  # cell G52
            length=1982,  # cell T52
            width_in=70,  # cell X52
            width_out=70,  # cell Y52
            height_in=55,  # cell Z52
            height_out=55,  # cell AA52
            mtop=3,  # cell AB52
            mbottom=3.5,  # cell AC52
            mleft=2.5,  # cell AD52
            mright=3,  # cell AE52
        ),
        "guide_06_02": guide(  # row 54
            x=30360.836074238643,  # cell H54
            y=97.97192314627146,  # cell I54
            z=-242.1691392655584,  # cell J54
            chi=-0.8242355981578852,  # cell F54
            phi=-0.41515089440691644,  # cell G54
            length=2000,  # cell T54
            width_in=70,  # cell X54
            width_out=70,  # cell Y54
            height_in=55,  # cell Z54
            height_out=55,  # cell AA54
            mtop=3,  # cell AB54
            mbottom=3.5,  # cell AC54
            mleft=2.5,  # cell AD54
            mright=3,  # cell AE54
        ),
        "guide_06_03": guide(  # row 56
            x=32362.57639085077,  # cell H56
            y=112.47828820342264,  # cell I56
            z=-270.96691605723163,  # cell J56
            chi=-0.7525442540421409,  # cell F56
            phi=-0.4438274320532142,  # cell G56
            length=2000,  # cell T56
            width_in=70,  # cell X56
            width_out=70,  # cell Y56
            height_in=55,  # cell Z56
            height_out=55,  # cell AA56
            mtop=3,  # cell AB56
            mbottom=3.5,  # cell AC56
            mleft=2.5,  # cell AD56
            mright=3,  # cell AE56
        ),
        "guide_06_04": guide(  # row 58
            x=34363.34377240512,  # cell H58
            y=127.97862955733744,  # cell I58
            z=-297.2474195217661,  # cell J58
            chi=-0.6808887197885923,  # cell F58
            phi=-0.4724896457546337,  # cell G58
            length=2000,  # cell T58
            width_in=70,  # cell X58
            width_out=70,  # cell Y58
            height_in=55,  # cell Z58
            height_out=55,  # cell AA58
            mtop=3,  # cell AB58
            mbottom=3.5,  # cell AC58
            mleft=2.5,  # cell AD58
            mright=3,  # cell AE58
        ),
        "guide_06_05": guide(  # row 60
            x=36365.13435204425,  # cell H60
            y=144.4884359197429,  # cell I60
            z=-321.0368758784358,  # cell J60
            chi=-0.609197375672848,  # cell F60
            phi=-0.5011661834009314,  # cell G60
            length=2000,  # cell T60
            width_in=70,  # cell X60
            width_out=70,  # cell Y60
            height_in=55,  # cell Z60
            height_out=55,  # cell AA60
            mtop=3,  # cell AB60
            mbottom=3.5,  # cell AC60
            mleft=2.5,  # cell AD60
            mright=3,  # cell AE60
        ),
        "guide_06_06": guide(  # row 62
            x=38365.94470731319,  # cell H62
            y=161.99120984653027,  # cell I62
            z=-342.31148223293405,  # cell J62
            chi=-0.5375418414192995,  # cell F62
            phi=-0.5298283971023509,  # cell G62
            length=2000,  # cell T62
            width_in=70,  # cell X62
            width_out=70,  # cell Y62
            height_in=55,  # cell Z62
            height_out=55,  # cell AA62
            mtop=3,  # cell AB62
            mbottom=3.5,  # cell AC62
            mleft=2.5,  # cell AD62
            mright=3,  # cell AE62
        ),
        "guide_06_07": guide(  # row 64
            x=40367.77101425058,  # cell H64
            y=180.50444097961105,  # cell I64
            z=-361.0924692474658,  # cell J64
            chi=-0.4658504973035552,  # cell F64
            phi=-0.5585049347486487,  # cell G64
            length=2000,  # cell T64
            width_in=70,  # cell X64
            width_out=70,  # cell Y64
            height_in=55,  # cell Z64
            height_out=55,  # cell AA64
            mtop=3,  # cell AB64
            mbottom=3.5,  # cell AC64
            mleft=2.5,  # cell AD64
            mright=3,  # cell AE64
        ),
        "guide_06_08": guide(  # row 66
            x=42368.60981467319,  # cell H66
            y=200.0096299502345,  # cell I66
            z=-377.36104532618344,  # cell J66
            chi=-0.3941949630500067,  # cell F66
            phi=-0.5871671484500681,  # cell G66
            length=2000,  # cell T66
            width_in=70,  # cell X66
            width_out=70,  # cell Y66
            height_in=55,  # cell Z66
            height_out=55,  # cell AA66
            mtop=3,  # cell AB66
            mbottom=3.5,  # cell AC66
            mleft=2.5,  # cell AD66
            mright=3,  # cell AE66
        ),
        "guide_06_09": guide(  # row 68
            x=44370.457313033396,  # cell H68
            y=220.52626731298074,  # cell I68
            z=-391.133445441332,  # cell J68
            chi=-0.3225036189342624,  # cell F68
            phi=-0.615843686096366,  # cell G68
            length=2000,  # cell T68
            width_in=70,  # cell X68
            width_out=70,  # cell Y68
            height_in=55,  # cell Z68
            height_out=55,  # cell AA68
            mtop=3,  # cell AB68
            mbottom=3.5,  # cell AC68
            mleft=2.5,  # cell AD68
            mright=3,  # cell AE68
        ),
        "guide_06_10": guide(  # row 70
            x=46371.31002997996,  # cell H70
            y=242.03385179298422,  # cell I70
            z=-402.39588941314656,  # cell J70
            chi=-0.25084808468071385,  # cell F70
            phi=-0.6445058997977854,  # cell G70
            length=2000,  # cell T70
            width_in=70,  # cell X70
            width_out=70,  # cell Y70
            height_in=55,  # cell Z70
            height_out=55,  # cell AA70
            mtop=3,  # cell AB70
            mbottom=3.5,  # cell AC70
            mleft=2.5,  # cell AD70
            mright=3,  # cell AE70
        ),
        "guide_06_11": guide(  # row 72
            x=48373.16418389658,  # cell H72
            y=264.55387483797375,  # cell I72
            z=-411.15961642229144,  # cell J72
            chi=-0.17915674056496955,  # cell F72
            phi=-0.6731824374440832,  # cell G72
            length=2000,  # cell T72
            width_in=70,  # cell X72
            width_out=70,  # cell Y72
            height_in=55,  # cell Z72
            height_out=55,  # cell AA72
            mtop=3,  # cell AB72
            mbottom=3.5,  # cell AC72
            mleft=2.5,  # cell AD72
            mright=3,  # cell AE72
        ),
        "guide_06_12": guide(  # row 74
            x=50374.01628882427,  # cell H74
            y=288.06383328750104,  # cell I74
            z=-417.4158577913398,  # cell J74
            chi=-0.10750120631142096,  # cell F74
            phi=-0.7018446511455027,  # cell G74
            length=2000,  # cell T74
            width_in=70,  # cell X74
            width_out=70,  # cell Y74
            height_in=55,  # cell Z74
            height_out=55,  # cell AA74
            mtop=3,  # cell AB74
            mbottom=3.5,  # cell AC74
            mleft=2.5,  # cell AD74
            mright=3,  # cell AE74
        ),
        "guide_06_13": guide(  # row 76
            x=52375.86256259569,  # cell H76
            y=312.58721946092004,  # cell I76
            z=-421.1708568390251,  # cell J76
            chi=-0.035809862195676695,  # cell F76
            phi=-0.7305211887918005,  # cell G76
            length=2000,  # cell T76
            width_in=70,  # cell X76
            width_out=70,  # cell Y76
            height_in=55,  # cell Z76
            height_out=55,  # cell AA76
            mtop=3,  # cell AB76
            mbottom=3.5,  # cell AC76
            mleft=2.5,  # cell AD76
            mright=3,  # cell AE76
        ),
    }
    return [to_bwc1, to_bwc2, to_end_of_bender]
