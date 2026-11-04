"""The straight guides of T-Rex, from guide 30 to the sample position.

Dimensions are from the Table of Optics V3.2 (the .xlsx file in ``reference/``), whose rows are
named in the comments. The values are the sheet's, not rounded. The m-values are those at the
entrance, which are the same as at the exit.

The axis is straight from ``arm_end_of_curve`` of :mod:`.curve` on, so nothing here is
rotated.

Guides 34 to 50a approximate the horizontal and vertical ellipses of the sheet with straight
tapered sections.
"""

from mcdgs.core.components import Arm, Component, Monitor_nD, Slit

from .curve import GUIDE_29_OUT_Z
from .sections import chain, gap, guide, split_after

ARM_MC_Z = 161.985  # m, cell B140
MC_GAP_LENGTH = 0.095  # m, cell D140
ARM_PSC_Z = 107.9  # m, cell B115
PSC_GAP_LENGTH = 0.2  # m, cell D115


def guides_straight() -> list[dict[str, Component]]:
    """Create the straight guides.

    Returns
    -------
    :
        Components by name, in beam order, in three pieces: up to ``arm_psc``, up to
        ``arm_mc``, and the rest.
    """
    guide_30_in_z = 96.00255  # m, cell B108

    sections = {
        "guide_30": guide(  # row 108
            z=guide_30_in_z,
            w1=0.06,  # cell E108
            h1=0.085,  # cell F108
            w2=0.06,  # cell G108
            h2=0.085,  # cell H108
            l=3.997,  # cell D108
            mleft=1,  # cell O108
            mright=1,  # cell P108
            mtop=2.5,  # cell Q108
            mbottom=2.5,  # cell R108
            at=[0, 0, guide_30_in_z - GUIDE_29_OUT_Z],
            relative="arm_end_of_curve",
        ),
        "guide_31": guide(  # row 110
            z=100.00255,  # cell B110
            w1=0.06,  # cell E110
            h1=0.085,  # cell F110
            w2=0.06,  # cell G110
            h2=0.085,  # cell H110
            l=3.997,  # cell D110
            mleft=1,  # cell O110
            mright=1,  # cell P110
            mtop=2.5,  # cell Q110
            mbottom=2.5,  # cell R110
        ),
        "guide_32": guide(  # row 112
            z=104.00255,  # cell B112
            w1=0.06,  # cell E112
            h1=0.085,  # cell F112
            w2=0.06,  # cell G112
            h2=0.085,  # cell H112
            l=3.89375,  # cell D112
            mleft=1,  # cell O112
            mright=1,  # cell P112
            mtop=2.5,  # cell Q112
            mbottom=2.5,  # cell R112
        ),
        "arm_psc": gap(  # row 115, the P-Chopper
            z=ARM_PSC_Z,
            length=PSC_GAP_LENGTH,
        ),
        "guide_33": guide(  # row 118
            z=108.1037,  # cell B118
            w1=0.06,  # cell E118
            h1=0.085,  # cell F118
            w2=0.06,  # cell G118
            h2=0.085,  # cell H118
            l=2.9943,  # cell D118
            mleft=1,  # cell O118
            mright=1,  # cell P118
            mtop=2.5,  # cell Q118
            mbottom=2.5,  # cell R118
        ),
        # row 119, the gap of beam monitor 2, overlaps guide 33 in the sheet and is left out
        "guide_34": guide(  # row 120
            z=111.098,  # cell B120
            w1=0.06012944068706193,  # cell E120
            h1=0.085,  # cell F120
            w2=0.0643426131,  # cell G120
            h2=0.085,  # cell H120
            l=4,  # cell D120
            mleft=1,  # cell O120
            mright=1,  # cell P120
            mtop=2.5,  # cell Q120
            mbottom=2.5,  # cell R120
        ),
        "guide_35": guide(  # row 121
            z=115.098,  # cell B121
            w1=0.0643426131,  # cell E121
            h1=0.085,  # cell F121
            w2=0.067328847516735621,  # cell G121
            h2=0.085,  # cell H121
            l=4,  # cell D121
            mleft=1,  # cell O121
            mright=1,  # cell P121
            mtop=2.5,  # cell Q121
            mbottom=2.5,  # cell R121
        ),
        "guide_36": guide(  # row 122
            z=119.098,  # cell B122
            w1=0.067328847516735621,  # cell E122
            h1=0.085,  # cell F122
            w2=0.069247059022427507,  # cell G122
            h2=0.085,  # cell H122
            l=4,  # cell D122
            mleft=1,  # cell O122
            mright=1,  # cell P122
            mtop=2.5,  # cell Q122
            mbottom=2.5,  # cell R122
        ),
        "guide_37": guide(  # row 123
            z=123.098,  # cell B123
            w1=0.069247059022427507,  # cell E123
            h1=0.085,  # cell F123
            w2=0.069542078013643575,  # cell G123
            h2=0.085,  # cell H123
            l=0.902,  # cell D123
            mleft=1,  # cell O123
            mright=1,  # cell P123
            mtop=2.5,  # cell Q123
            mbottom=2.5,  # cell R123
        ),
        "guide_38": guide(  # row 124
            z=124,  # cell B124
            w1=0.069542078013643575,  # cell E124
            h1=0.084998466775413007,  # cell F124
            w2=0.070265720819484088,  # cell G124
            h2=0.084562370357825117,  # cell H124
            l=4,  # cell D124
            mleft=1,  # cell O124
            mright=1,  # cell P124
            mtop=2.5,  # cell Q124
            mbottom=2.5,  # cell R124
        ),
        "guide_39": guide(  # row 125
            z=128,  # cell B125
            w1=0.070265720819484088,  # cell E125
            h1=0.084562370357825117,  # cell F125
            w2=0.070051567125885967,  # cell G125
            h2=0.083339719228873051,  # cell H125
            l=4,  # cell D125
            mleft=1,  # cell O125
            mright=1,  # cell P125
            mtop=2.5,  # cell Q125
            mbottom=2.5,  # cell R125
        ),
        "guide_40": guide(  # row 126
            z=132,  # cell B126
            w1=0.070051567125885967,  # cell E126
            h1=0.083339719228873051,  # cell F126
            w2=0.0688908718,  # cell G126
            h2=0.081295032542039323,  # cell H126
            l=4,  # cell D126
            mleft=1,  # cell O126
            mright=1,  # cell P126
            mtop=2.5,  # cell Q126
            mbottom=2.5,  # cell R126
        ),
        "guide_41": guide(  # row 127
            z=136,  # cell B127
            w1=0.0688908718,  # cell E127
            h1=0.081295032542039323,  # cell F127
            w2=0.0667342641,  # cell G127
            h2=0.078363990615112868,  # cell H127
            l=4,  # cell D127
            mleft=1,  # cell O127
            mright=1,  # cell P127
            mtop=2.5,  # cell Q127
            mbottom=2.5,  # cell R127
        ),
        "guide_42": guide(  # row 128
            z=140,  # cell B128
            w1=0.0667342641,  # cell E128
            h1=0.078363990615112868,  # cell F128
            w2=0.0634803231,  # cell G128
            h2=0.0744419702,  # cell H128
            l=4,  # cell D128
            mleft=1,  # cell O128
            mright=1,  # cell P128
            mtop=3,  # cell Q128
            mbottom=3,  # cell R128
        ),
        "guide_43": guide(  # row 129
            z=144,  # cell B129
            w1=0.0634803231,  # cell E129
            h1=0.0744419702,  # cell F129
            w2=0.058947607768968737,  # cell G129
            h2=0.069361069960787788,  # cell H129
            l=4,  # cell D129
            mleft=1.5,  # cell O129
            mright=1.5,  # cell P129
            mtop=3,  # cell Q129
            mbottom=3,  # cell R129
        ),
        "guide_44": guide(  # row 130
            z=148,  # cell B130
            w1=0.058947607768968737,  # cell E130
            h1=0.069361069960787788,  # cell F130
            w2=0.053665251983516185,  # cell G130
            h2=0.063732637664680929,  # cell H130
            l=3.509,  # cell D130
            mleft=1.5,  # cell O130
            mright=1.5,  # cell P130
            mtop=3,  # cell Q130
            mbottom=3,  # cell R130
        ),
        "guide_45": guide(  # row 132
            z=151.596,  # cell B132
            w1=0.053665251983516185,  # cell E132
            h1=0.063732637664680929,  # cell F132
            w2=0.045385122075640119,  # cell G132
            h2=0.055334681646321194,  # cell H132
            l=4,  # cell D132
            mleft=1.5,  # cell O132
            mright=1.5,  # cell P132
            mtop=3,  # cell Q132
            mbottom=3,  # cell R132
        ),
        # the sheet has no guide 46
        "guide_47": guide(  # row 133
            z=155.596,  # cell B133
            w1=0.045385122075640119,  # cell E133
            h1=0.055334681646321194,  # cell F133
            w2=0.038917728211423487,  # cell G133
            h2=0.049103212537757955,  # cell H133
            l=2.399,  # cell D133
            mleft=2,  # cell O133
            mright=2,  # cell P133
            mtop=4,  # cell Q133
            mbottom=4,  # cell R133
        ),
        # guides 48a and 48c, the polarizers, are not modeled: 48b is the plain guide in their slot
        "guide_48b": guide(  # row 136
            z=158,  # cell B136
            w1=0.038917728211423487,  # cell E136
            h1=0.049103212537757955,  # cell F136
            w2=0.035667285784256798,  # cell G136
            h2=0.046091729354631667,  # cell H136
            l=1,  # cell D136
            mleft=2,  # cell O136
            mright=2,  # cell P136
            mtop=4.5,  # cell Q136
            mbottom=4.5,  # cell R136
        ),
        "guide_49": guide(  # row 139
            z=159.003,  # cell B139
            w1=0.035656965524720398,  # cell E139
            h1=0.04608231057651308,  # cell F139
            w2=0.022433235171385018,  # cell G139
            h2=0.034986956776210739,  # cell H139
            l=2.982,  # cell D139
            mleft=2,  # cell O139
            mright=2,  # cell P139
            mtop=4.5,  # cell Q139
            mbottom=4.5,  # cell R139
        ),
        "arm_mc": gap(  # row 140, the M-Chopper
            z=ARM_MC_Z,
            length=MC_GAP_LENGTH,
        ),
        "guide_50a": guide(  # row 142
            z=162.083,  # cell B142
            w1=0.022433235171385018,  # cell E142
            h1=0.034986956776210739,  # cell F142
            w2=0.014184368544301893,  # cell G142
            h2=0.029555541228075157,  # cell H142
            l=1.003,  # cell D142
            mleft=3.5,  # cell O142
            mright=3.5,  # cell P142
            mtop=4.5,  # cell Q142
            mbottom=4.5,  # cell R142
        ),
    }
    return split_after(chain(sections), "arm_psc", "arm_mc")


def slit_mc() -> Slit:
    """Create the slit after the last guide.

    Returns
    -------
    :
        A new slit.
    """
    slit_mc_z = 163.105  # m, cell AB145

    return Slit(
        xwidth=0.0205,  # FIXME, not in Table of Optics v3.2
        yheight=0.0345,  # FIXME, not in Table of Optics v3.2
        at=[0, 0, slit_mc_z - ARM_MC_Z],
        relative="arm_mc",
    )


def sample_pos() -> Arm:
    """Create the arm at the sample position.

    Returns
    -------
    :
        A new arm.
    """
    sample_z = 163.8  # m, cell C148

    return Arm(at=[0, 0, sample_z - ARM_MC_Z], relative="arm_mc")


def sample_monitor() -> Monitor_nD:
    """Create the wavelength monitor at the sample position.

    It records the beam reaching the sample and leaves the neutrons unchanged. It is at
    ``sample_pos``, before the sample.

    Returns
    -------
    :
        A new monitor.
    """
    return Monitor_nD(
        xwidth=0.1,  # FIXME, not in Table of Optics v3.2: wider than the beam
        yheight=0.1,  # FIXME, not in Table of Optics v3.2: wider than the beam
        restore_neutron=1,
        options="lambda limits=[0.5 5.3] bins=12",  # FIXME, not in Table of Optics v3.2
        filename="sample_monitor",
        relative="sample_pos",
    )
