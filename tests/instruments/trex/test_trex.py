import math
from itertools import pairwise

import numpy as np
import pytest

from mcdgs.core.components import Guide_gravity
from mcdgs.core.geometry import place_absolute
from mcdgs.instruments.trex.guides.extraction import half_taper_angle
from mcdgs.instruments.trex.model import TRex
from mcdgs.samples import aluminum


def test_the_components_are_in_beam_order():
    assert list(TRex().components) == [
        "origin",
        "source",
        "arm_NBOA",
        "guide_1",
        "guide_2",
        "guide_3",
        "arm_guide_4",
        "guide_4a",
        "bender",
        "guide_4b",
        "arm_start_of_curve",
        "guide_5",
        "guide_6",
        "guide_7",
        "guide_8",
        "guide_9",
        "guide_10",
        "guide_11",
        "guide_12",
        "guide_13a",
        "arm_bwc1",
        "guide_13b",
        "guide_14a",
        "guide_14b",
        "guide_15a",
        "arm_bwc2",
        "guide_15b",
        "arm_bm1",
        *[f"guide_{number}" for number in range(16, 30)],
        "arm_end_of_curve",
        *[f"guide_{number}" for number in range(30, 33)],
        "arm_psc",
        *[f"guide_{number}" for number in range(33, 46)],
        "guide_47",
        "guide_48b",
        "guide_49",
        "arm_mc",
        "guide_50a",
        "slit_mc",
        "sample_pos",
    ]


def test_the_source_is_that_of_trex_primary():
    source = TRex().source
    assert (source.sector, source.beamline, source.acc_power) == ("W", 7, 2)
    assert (source.yheight, source.cold_frac, source.dist) == (0.03, 0.3, 2)
    assert (source.focus_xw, source.focus_yh) == (0.095, 0.035)
    assert (source.Lmin, source.Lmax, source.n_pulses) == (0.1, 10.0, 1)


def test_the_sample_position_is_on_the_axis_after_the_m_chopper():
    trex = TRex()
    assert trex.sample_pos.relative == "arm_mc"
    assert trex.sample_pos.at == pytest.approx([0, 0, 163.8 - 161.985], abs=1e-12)
    assert trex.sample_pos.rotated == [0, 0, 0]


def test_without_a_sample_there_is_none():
    trex = TRex()
    assert trex.sample is None
    assert "sample" not in trex.components


def test_a_sample_goes_after_the_sample_position():
    trex = TRex(sample=aluminum())
    assert list(trex.components)[-2:] == ["sample_pos", "sample"]
    assert trex.sample.relative == "sample_pos"
    assert trex.positions()["sample"] == trex.positions()["sample_pos"]


def test_a_sample_keeps_its_own_placement():
    trex = TRex(sample=aluminum(relative="guide_1"))
    assert trex.sample.relative == "guide_1"


def test_the_given_sample_is_not_changed():
    crystal = aluminum()
    trex = TRex(sample=crystal)
    assert trex.sample is not crystal
    assert (crystal.name, crystal.relative) == (None, None)


def test_the_components_are_those_of_the_table_of_optics():
    trex = TRex()
    assert trex.arm_NBOA.at == pytest.approx([-0.002022816656807591, 0, 1.895927], abs=1e-12)
    assert trex.arm_NBOA.rotated == pytest.approx([0, -0.21499811064881996, 0], abs=1e-12)
    assert trex.arm_NBOA.relative is None
    assert trex.guide_1 == Guide_gravity(
        w1=0.089955,
        h1=0.029992,
        w2=0.082455,
        h2=0.034937,
        l=0.9994,
        mleft=1.5,
        mright=1.5,
        mtop=4,
        mbottom=4,
        at=[0, 0, 0],
        rotated=[0, 0, 0],
        relative="arm_NBOA",
    )


def test_the_guides_start_at_the_lab_z_of_the_sheet():
    positions = TRex().positions()
    assert positions["guide_1"][2] == pytest.approx(1.895927, abs=1e-9)
    assert positions["guide_2"][2] == pytest.approx(2.896327, abs=1e-9)
    assert positions["guide_3"][2] == pytest.approx(3.896727, abs=1e-9)
    assert positions["guide_4a"][2] == pytest.approx(5.40065, abs=1e-9)
    assert positions["bender"][2] == pytest.approx(5.40065, abs=1e-9)
    assert positions["guide_4b"][2] == pytest.approx(5.45065, abs=1e-9)


def test_the_guides_hang_off_the_arm():
    trex = TRex()
    assert trex.guide_2.relative == trex.guide_3.relative == "arm_NBOA"
    assert trex.guide_2.at[:2] == trex.guide_3.at[:2] == [0, 0]


def test_guide_3_keeps_the_m_of_mctrex_where_the_sheet_says_issue():
    guide_3 = TRex().guide_3
    assert (guide_3.mtop, guide_3.mbottom) == (4, 4)


def test_the_guide_starts_on_the_bisector_of_the_fixed_wall():
    trex = TRex()
    x, y, z = trex.positions()["guide_1"]
    assert x == pytest.approx(-0.0020228, abs=1e-6)
    assert (y, z) == (0.0, 1.895927)
    assert x - trex.guide_1.w1 / 2 == pytest.approx(-0.047, abs=1e-6)  # the -x wall


def test_the_tilt_is_half_the_taper_of_the_sheet():
    # cell N13 is atan(dw / L): half of it differs from atan(dw / L / 2) by about 3e-6 deg
    assert TRex().arm_NBOA.rotated[1] == pytest.approx(-0.42999016682346453 / 2, abs=1e-5)


def test_the_centre_line_tilts_towards_minus_x():
    trex = TRex()
    instr_param_dict = {p.name: p.value for p in trex.instr_params}
    (position, rot) = place_absolute(trex.components, instr_param_dict)["guide_1"]
    axis = rot[2]  # the guide's z axis in the lab frame
    at_exit = position + axis / axis[2] * (5.375927 - 1.895927)
    assert at_exit[0] == pytest.approx(-0.047 + 0.063838 / 2, abs=1e-6)
    assert at_exit[1] == pytest.approx(0)


def test_guide_4_tilts_by_half_its_own_taper():
    # cell N25, with the same difference to atan(dw / L / 2) as for cell N13
    trex = TRex()
    tilt = trex.arm_guide_4.rotated[1]
    assert tilt == pytest.approx(-0.42669474420770109 / 2, abs=1e-5)
    guide_4b = trex.guide_4b
    assert half_taper_angle(guide_4b.w1, guide_4b.w2, guide_4b.l) == pytest.approx(tilt, abs=1e-12)


def test_guides_4a_4b_and_the_bender_hang_off_one_arm():
    trex = TRex()
    assert trex.guide_4a.relative == trex.bender.relative == trex.guide_4b.relative == "arm_guide_4"
    assert trex.guide_4a.at == trex.bender.at == [0, 0, 0]
    assert trex.guide_4b.at[:2] == [0, 0]
    assert trex.arm_guide_4.when is None  # guide 4a and the bender need it in either mode


@pytest.mark.parametrize("name", ["guide_4a", "guide_4b"])
def test_the_minus_x_wall_of_guide_4_is_the_fixed_wall(name):
    trex = TRex()
    instr_param_dict = {p.name: p.value for p in trex.instr_params}
    (position, rot) = place_absolute(trex.components, instr_param_dict)[name]
    guide = getattr(trex, name)
    # the rows of rot are the guide's axes in the lab frame
    at_entrance = position - rot[0] * guide.w1 / 2
    at_exit = position - rot[0] * guide.w2 / 2 + rot[2] * guide.l
    assert at_entrance[0] == pytest.approx(-0.047, abs=1e-8)
    assert at_exit[0] == pytest.approx(-0.047, abs=1e-8)


def test_the_bender_selects_guide_4a_or_guide_4b():
    trex = TRex()
    assert trex.params.bender_in is False
    assert trex.guide_4a.when == "bender_in==0"
    assert trex.bender.when == trex.guide_4b.when == "bender_in==1"


def test_every_when_names_an_instrument_parameter():
    trex = TRex()
    names = {p.name for p in trex.instr_params}
    whens = [c.when for c in trex.components.values() if c.when is not None]
    assert whens
    assert all(when.split("==")[0] in names for when in whens)


def test_guides_4a_and_4b_share_their_exit():
    trex = TRex()
    assert (trex.guide_4a.w2, trex.guide_4a.h2) == (trex.guide_4b.w2, trex.guide_4b.h2)
    guide_4a_exit = trex.positions()["guide_4a"][2] + trex.guide_4a.l
    guide_4b_exit = trex.positions()["guide_4b"][2] + trex.guide_4b.l
    assert guide_4a_exit == pytest.approx(guide_4b_exit, abs=1e-9)


def test_the_bender_is_rotated_by_its_full_deflection():
    trex = TRex()
    assert trex.bender.rotated == [0, "b_rot", 0]
    assert trex.params.b_rot == pytest.approx(-0.397881, abs=1e-6)
    assert {p.name: p.value for p in trex.instr_params} == {
        "bender_in": 0,
        "b_rot": trex.params.b_rot,
    }


def test_every_mirror_of_the_bender_is_the_same():
    bender = TRex().bender
    mirrors = [
        getattr(bender, f"r{wall}{spin}Par")
        for wall in ("Top", "Bot", "Left", "Right")
        for spin in ("Up", "Down")
    ]
    assert all(mirror == [0.99, 0.0219, 3.02, 4.0, 0.003] for mirror in mirrors)
    assert mirrors[0] is not mirrors[1]


def lab_placement(trex: TRex) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    instr_param_dict = {p.name: p.value for p in trex.instr_params}
    return place_absolute(trex.components, instr_param_dict)


def bearing(rot: np.ndarray) -> float:
    """Angle of a component's z axis to the lab z axis, in deg. The rows of rot are its axes."""
    return math.degrees(math.atan2(rot[2][0], rot[2][2]))


def test_the_curve_starts_at_the_lab_z_of_the_sheet_with_the_tilt_taken_out():
    trex = TRex()
    assert trex.arm_start_of_curve.relative == "arm_guide_4"
    position, rot = lab_placement(trex)["arm_start_of_curve"]
    assert position[2] == pytest.approx(5.916, abs=1e-9)
    assert bearing(rot) == pytest.approx(0, abs=1e-12)


def test_guide_5_turns_by_half_its_length_from_the_start_of_the_curve():
    trex = TRex()
    assert trex.guide_5.relative == "arm_start_of_curve"
    assert trex.guide_5.at == [0, 0, 0]
    assert bearing(lab_placement(trex)["guide_5"][1]) == pytest.approx(
        math.degrees(2.7625 / 2 / 12000), abs=1e-9
    )


def test_the_guides_are_those_of_the_table_of_optics():
    trex = TRex()
    assert trex.guide_11 == Guide_gravity(  # row 60, the first with different side m-values
        w1=0.06,
        h1=0.078035464780223663,
        w2=0.06,
        h2=0.080841682168482322,
        l=3.6349,
        mleft=1.5,
        mright=2.5,
        mtop=2.5,
        mbottom=2.5,
        at=trex.guide_11.at,
        rotated=trex.guide_11.rotated,
        relative="guide_10",
    )
    assert trex.guide_11.at == pytest.approx([0, 0, 24.3017 - 22.7007], abs=1e-12)
    assert trex.guide_11.rotated == pytest.approx(
        [0, math.degrees(((1.588 + 3.6349) / 2 + 0.013) / 12000), 0], abs=1e-12
    )
    assert trex.guide_50a == Guide_gravity(  # row 142
        w1=0.022433235171385018,
        h1=0.034986956776210739,
        w2=0.014184368544301893,
        h2=0.029555541228075157,
        l=1.003,
        mleft=3.5,
        mright=3.5,
        mtop=4.5,
        mbottom=4.5,
        at=trex.guide_50a.at,
        rotated=[0, 0, 0],
        relative="arm_mc",
    )
    assert trex.guide_50a.at == pytest.approx([0, 0, 162.083 - 161.985], abs=1e-12)


def test_each_guide_follows_the_one_before_on_the_axis():
    trex = TRex()
    components = list(trex.components)
    start = components.index("guide_6")
    stop = components.index("slit_mc")
    for previous, name in pairwise(components[start - 1 : stop]):
        assert trex.components[name].relative == previous
        assert trex.components[name].at[:2] == [0, 0]
        assert trex.components[name].at[2] >= 0


@pytest.mark.parametrize(
    ("name", "following", "distance"),
    [
        ("arm_bwc1", "guide_13b", 0.03),  # cell B66 minus cell B65
        ("arm_bwc2", "guide_15b", 0.03),  # cell B74 minus cell B73
        ("arm_bm1", "guide_16", 0.02369),  # cell B80 minus cell B77
        ("arm_psc", "guide_33", 0.2037),  # cell B118 minus cell B115
        ("arm_mc", "guide_50a", 0.098),  # cell B142 minus cell B140
    ],
)
def test_the_next_guide_is_on_the_arm_of_the_gap(name, following, distance):
    guide = TRex().components[following]
    assert guide.relative == name
    assert guide.at[2] == pytest.approx(distance, abs=1e-9)


def test_the_length_along_the_axis_is_that_of_the_sheet():
    trex = TRex()
    positions = lab_placement(trex)
    names = list(trex.components)
    path = sum(
        float(np.linalg.norm(positions[name][0] - positions[previous][0]))
        for previous, name in pairwise(
            names[names.index("guide_5") : names.index("sample_pos") + 1]
        )
    )
    assert path == pytest.approx(163.8 - 5.916, abs=1e-9)


def test_the_axis_leaves_the_curve_on_its_tangent():
    trex = TRex()
    placement = lab_placement(trex)
    end_of_curve = 95.99955 - 5.916  # cell C106 minus cell B32
    assert bearing(placement["arm_end_of_curve"][1]) == pytest.approx(
        math.degrees(end_of_curve / 12000), abs=1e-9
    )
    assert trex.guide_30.relative == "arm_end_of_curve"
    assert trex.guide_30.at == pytest.approx([0, 0, 0.003], abs=1e-9)  # cell D107
    assert trex.guide_30.rotated == [0, 0, 0]
    assert bearing(placement["sample_pos"][1]) == pytest.approx(
        bearing(placement["guide_30"][1]), abs=1e-9
    )


def test_each_guide_on_the_curve_lies_at_the_angle_of_the_arc_at_its_midpoint():
    trex = TRex()
    placement = lab_placement(trex)
    midpoint = 37.94924 + 1.9908 / 2  # guide 15a, row 72
    assert bearing(placement["guide_15a"][1]) == pytest.approx(
        math.degrees((midpoint - 5.916) / 12000), abs=1e-9
    )


def test_the_sample_position_agrees_with_the_design_model():
    # the estimate from trex_guide_anyshape.ipynb of mctrex, its z is lab z and not on the axis
    position, rot = lab_placement(TRex())["sample_pos"]
    offset = rot @ (np.array([0.8301227276790154, 0, 163.8]) - position)
    assert offset[0] == pytest.approx(0, abs=1e-4)
    assert bearing(rot) == pytest.approx(0.430122381481144, abs=1e-5)


def test_the_slit_is_after_guide_50a():
    trex = TRex()
    assert trex.slit_mc.relative == "arm_mc"
    assert trex.slit_mc.at[2] == pytest.approx(163.105 - 161.985, abs=1e-12)
    assert trex.slit_mc.at[2] > trex.guide_50a.at[2] + trex.guide_50a.l


def test_two_instruments_share_no_component():
    first, second = TRex(), TRex()
    assert first == second
    assert first.guide_1 is not second.guide_1
    first.guide_1.w1 = 0.1
    assert second.guide_1.w1 == pytest.approx(0.089955)
