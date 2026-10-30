import pytest

from mcdgs.core.components import Guide_gravity
from mcdgs.core.geometry import place_absolute
from mcdgs.instruments.trex.guides_extraction import (
    GUIDE_1_IN_Z,
    GUIDE_3_OUT_WIDTH,
    GUIDE_3_OUT_Z,
    NBOA_TAPER_ANGLE,
)
from mcdgs.instruments.trex.trex import TRex
from mcdgs.samples import aluminum


def test_the_components_are_in_beam_order():
    assert list(TRex().components) == [
        "origin",
        "source",
        "arm_NBOA",
        "guide_1",
        "guide_2",
        "guide_3",
        "guide_4a",
        "bender",
        "guide_4b",
        "sample_pos",
    ]


def test_the_source_is_that_of_trex_primary():
    source = TRex().source
    assert (source.sector, source.beamline, source.acc_power) == ("W", 7, 2)
    assert (source.yheight, source.cold_frac, source.dist) == (0.03, 0.3, 2)
    assert (source.focus_xw, source.focus_yh) == (0.095, 0.035)
    assert (source.Lmin, source.Lmax, source.n_pulses) == (0.1, 10.0, 1)


def test_the_sample_position_is_absolute():
    trex = TRex()
    assert trex.sample_pos.relative is None
    assert trex.sample_pos.rotated == [0, 0.430122381481144, 0]
    assert trex.positions()["sample_pos"] == pytest.approx([0.830122728, 0, 163.8], abs=1e-9)


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
    assert (y, z) == (0.0, GUIDE_1_IN_Z)
    assert x - trex.guide_1.w1 / 2 == pytest.approx(-0.047, abs=1e-6)  # the -x wall


def test_the_tilt_is_half_the_taper_of_the_sheet():
    # cell N13 is atan(dw / L): half of it differs from atan(dw / L / 2) by about 3e-6 deg
    assert NBOA_TAPER_ANGLE == pytest.approx(0.42999016682346453 / 2, abs=1e-5)


def test_the_centre_line_tilts_towards_minus_x():
    trex = TRex()
    instr_param_dict = {p.name: p.value for p in trex.instr_params}
    (position, rot) = place_absolute(trex.components, instr_param_dict)["guide_1"]
    axis = rot[2]  # the guide's z axis in the lab frame
    at_exit = position + axis / axis[2] * (GUIDE_3_OUT_Z - GUIDE_1_IN_Z)
    assert at_exit[0] == pytest.approx(-0.047 + GUIDE_3_OUT_WIDTH / 2, abs=1e-6)
    assert at_exit[1] == pytest.approx(0)


def test_the_bender_selects_guide_4a_or_guide_4b():
    trex = TRex()
    assert trex.params.bender == 0
    assert trex.guide_4a.when == "bender==0"
    assert trex.bender.when == trex.guide_4b.when == "bender==1"


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
        "bender": 0,
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


def test_two_instruments_share_no_component():
    first, second = TRex(), TRex()
    assert first == second
    assert first.guide_1 is not second.guide_1
    first.guide_1.w1 = 0.1
    assert second.guide_1.w1 == pytest.approx(0.089955)
