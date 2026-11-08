import math
from itertools import pairwise

import numpy as np
import pytest

from mcdgs.core.components import Arm, Guide_gravity, NXdisk_chopper
from mcdgs.core.geometry import place_absolute
from mcdgs.instruments.trex.choppers import (
    BEAM_ANGLE,
    DISK_DISTANCE,
    MODE_SLIT,
    SLIT_EDGES,
    chopper_timings,
    slit_angle,
)
from mcdgs.instruments.trex.guides.extraction import half_taper_angle
from mcdgs.instruments.trex.guides.sections import split_after
from mcdgs.instruments.trex.model import TRex
from mcdgs.samples import aluminum


def test_the_components_are_in_beam_order():
    assert list(TRex().components) == [
        "origin",
        "source",
        "arm_nboa",
        "guide_1",
        "guide_2",
        "guide_3",
        "arm_guide_4",
        "guide_4a",
        "bender",
        "guide_4b",
        "beam_monitor",
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
        "bw1",
        "guide_13b",
        "guide_14a",
        "guide_14b",
        "guide_15a",
        "arm_bwc2",
        "bw2",
        "guide_15b",
        "arm_bm1",
        *[f"guide_{number}" for number in range(16, 30)],
        "arm_end_of_curve",
        *[f"guide_{number}" for number in range(30, 33)],
        "arm_psc",
        "ps1",
        "ps2",
        *[f"guide_{number}" for number in range(33, 46)],
        "guide_47",
        "guide_48b",
        "guide_49",
        "arm_mc",
        "mc1",
        "mc2",
        "guide_50a",
        "slit_mc",
        "pre_sample_monitor",
        "beamstop_monitor",
        "sample_pos",
        "sample_monitor",
        "beamstop",
        "banana_detector",
    ]


# half the band: (h / m_n) (1 / 14 Hz + 3.5 ms) / 163.8 m, over 2
HALF_BAND = 0.5 * 3956.034006 * (1 / 14 + 3.5e-3) / 163.8  # Å


def test_the_source_is_that_of_trex_primary():
    source = TRex().source
    assert (source.sector, source.beamline, source.acc_power) == ("W", 7, 2)
    assert (source.yheight, source.cold_frac) == (0.03, 0.3)
    assert source.n_pulses == 1


def test_the_source_focuses_on_the_whole_entrance_of_guide_1():
    trex = TRex()
    source, guide_1, arm = trex.source, trex.guide_1, trex.components["arm_nboa"]
    x, _, z = arm.at
    assert source.dist == z  # the plane of the entrance
    assert source.focus_yh == guide_1.h1
    # centered on the axis, the rectangle reaches both walls of the entrance, off the axis
    half_width = guide_1.w1 * math.cos(math.radians(arm.rotated[1])) / 2
    walls = (x - half_width, x + half_width)
    assert walls == pytest.approx((-0.047, 0.042954), abs=1e-6)
    assert source.focus_xw / 2 == pytest.approx(max(abs(wall) for wall in walls))


def test_the_band_of_the_source_is_centered_on_lambda0():
    source = TRex().source
    assert (source.Lmin, source.Lmax) == pytest.approx((2.0 - HALF_BAND, 2.0 + HALF_BAND))


def test_a_new_lambda0_moves_the_band_of_the_source_and_the_time_and_wavelength_monitors():
    trex = TRex()
    names = list(trex.components)
    trex.update_params(lambda0=4.0)
    assert (trex.source.Lmin, trex.source.Lmax) == pytest.approx((4.0 - HALF_BAND, 4.0 + HALF_BAND))
    fresh = TRex(params=trex.params)
    for name in ("beam_monitor", "pre_sample_monitor", "beamstop_monitor"):
        assert getattr(trex, name) == getattr(fresh, name)
        assert getattr(trex, name).bin_limits["lambda"][0] == pytest.approx(4.0 - HALF_BAND)
    assert list(trex.components) == names  # in their place


def test_the_beam_monitor_is_in_the_middle_of_its_gap():
    position = TRex().positions()["beam_monitor"]
    assert position[2] == pytest.approx(5.892 + 0.02 / 2, abs=1e-9)  # cells B29 and D29


def test_the_sample_position_is_on_the_axis_after_the_m_chopper():
    trex = TRex()
    assert trex.sample_pos.relative == "arm_mc"
    assert trex.sample_pos.at == pytest.approx([0, 0, 163.8 - 161.985], abs=1e-12)
    assert trex.sample_pos.rotated == [0, 0, 0]


@pytest.mark.parametrize(
    ("name", "arm", "z"),
    [
        ("bw1", "arm_bwc1", 0.015),
        ("bw2", "arm_bwc2", 0.015),
        ("ps1", "arm_psc", 0.05),
        ("ps2", "arm_psc", 0.15),
        ("mc1", "arm_mc", 0.0425),
        ("mc2", "arm_mc", 0.0525),
    ],
)
def test_each_chopper_is_in_the_gap_of_its_arm(name, arm, z):
    chopper = getattr(TRex(), name)
    assert chopper.relative == arm
    assert chopper.at == pytest.approx([0, 0, z], abs=1e-12)


@pytest.mark.parametrize(
    ("name", "distance"),
    [  # the positions in column AB of the sheet, two disks around each P and M position
        ("bw1", 31.96404),  # row 65
        ("bw2", 39.95504),  # row 73
        ("ps1", 107.95),  # row 115
        ("ps2", 108.05),
        ("mc1", 162.0275),  # row 140
        ("mc2", 162.0375),
    ],
)
def test_each_chopper_is_at_its_distance_from_the_moderator(name, distance):
    assert DISK_DISTANCE[name] == pytest.approx(distance, abs=1e-9)


def test_the_choppers_turn_at_the_frequencies_of_trex_primary():
    timings = chopper_timings(lambda0=2.4, rrm=4, source_frequency=14, mode="FF", t0_offset=0.0017)
    assert {name: timing.nu for name, timing in timings.items()} == pytest.approx(
        {"bw1": 14, "bw2": 14, "ps1": 42, "ps2": -42, "mc1": 56, "mc2": 56}
    )


def test_in_rr_each_zero_mark_passes_the_beam_when_the_central_wavelength_arrives():
    timings = chopper_timings(lambda0=2.4, rrm=4, source_frequency=14, mode="RR", t0_offset=0.0017)
    for name, timing in timings.items():
        arrival = 0.0017 + DISK_DISTANCE[name] * 2.4 / 3956.034
        assert timing.delay == pytest.approx(arrival, rel=1e-6)


@pytest.mark.parametrize(
    ("mode", "angles"),
    [  # FF: the shifts of the slit edges in trex_primary.txt's HF mode
        ("FF", {"bw1": 0, "bw2": 0, "ps1": 55, "ps2": 55, "mc1": 175, "mc2": 175}),
        ("RR", {"bw1": 0, "bw2": 0, "ps1": 0, "ps2": 0, "mc1": 0, "mc2": 0}),
    ],
)
def test_the_slit_angles_are_the_shifts_of_trex_primary(mode, angles):
    assert {name: slit_angle(name, mode) for name in angles} == pytest.approx(angles)


@pytest.mark.parametrize("mode", ["FF", "RR"])
def test_the_slit_of_the_mode_is_centered_on_the_beam_when_the_central_wavelength_arrives(mode):
    timings = chopper_timings(lambda0=2.4, rrm=4, source_frequency=14, mode=mode, t0_offset=0.0017)
    for name, timing in timings.items():
        arrival = 0.0017 + DISK_DISTANCE[name] * 2.4 / 3956.034
        # as NXdisk_chopper: the disk angle on the beam is beam_angle minus the turn of the mark
        mark = 360 * timing.nu * (arrival - timing.delay)
        on_beam = BEAM_ANGLE[name] - mark
        edges = SLIT_EDGES[name]
        slit = MODE_SLIT[mode][name]
        center = (edges[2 * slit] + edges[2 * slit + 1]) / 2
        off_center = (on_beam - center + 180) % 360 - 180  # in [-180, 180): 359.9... is -0.0...
        assert off_center == pytest.approx(0, abs=1e-3), name


def test_the_mode_selects_the_wide_slit_for_flux_and_the_narrow_one_for_resolution():
    for name in ("ps1", "ps2", "mc1", "mc2"):
        edges = SLIT_EDGES[name]

        def width(mode, edges=edges, name=name):
            slit = MODE_SLIT[mode][name]
            return edges[2 * slit + 1] - edges[2 * slit]

        assert width("FF") > width("RR"), name


def test_an_unknown_mode_raises():
    with pytest.raises(ValueError, match="mode 'HF' is not 'FF' or 'RR'"):
        chopper_timings(lambda0=2.4, rrm=4, source_frequency=14, mode="HF", t0_offset=0.0017)


def test_each_chopper_takes_its_timing_from_the_instrument_parameters():
    trex = TRex()
    instr_param_dict = {p.name: p.value for p in trex.instr_params}
    timings = chopper_timings(
        trex.params.lambda0,
        trex.params.rrm,
        trex.source.frequency,
        trex.params.mode,
        trex.params.t0_offset,
    )
    for name, timing in timings.items():
        chopper = getattr(trex, name)
        assert (chopper.nu, chopper.delay) == (f"nu_{name}", f"delay_{name}")
        assert (instr_param_dict[chopper.nu], instr_param_dict[chopper.delay]) == timing


def test_the_timing_follows_the_control_parameters():
    trex = TRex()
    trex.update_params(lambda0=4.0, rrm=6)
    instr_param_dict = {p.name: p.value for p in trex.instr_params}
    assert instr_param_dict["nu_mc1"] == pytest.approx(84)
    assert instr_param_dict["delay_bw1"] == chopper_timings(4.0, 6, 14, "FF", 0.0017)["bw1"].delay


def test_the_t0_offset_delays_every_chopper_by_as_much():
    trex = TRex()
    before = {p.name: p.value for p in trex.instr_params}
    trex.update_params(t0_offset=trex.params.t0_offset + 0.001)
    after = {p.name: p.value for p in trex.instr_params}
    for name in ("bw1", "bw2", "ps1", "ps2", "mc1", "mc2"):
        assert after[f"delay_{name}"] - before[f"delay_{name}"] == pytest.approx(0.001)
        assert after[f"nu_{name}"] == before[f"nu_{name}"]


def test_the_default_mode_is_ff_and_rr_moves_only_the_pulse_shaping_and_monochromatic_choppers():
    trex = TRex()
    assert trex.params.mode == "FF"
    ff = {p.name: p.value for p in trex.instr_params}
    trex.update_params(mode="RR")
    rr = {p.name: p.value for p in trex.instr_params}
    changed = {name for name in ff if ff[name] != rr[name]}
    assert changed == {"delay_ps1", "delay_ps2", "delay_mc1", "delay_mc2"}


def test_a_cut_keeps_every_component_in_order():
    components = {name: Arm() for name in ("a", "b", "c", "d")}
    pieces = split_after(components, "a", "c")
    assert [list(piece) for piece in pieces] == [["a"], ["b", "c"], ["d"]]


def test_without_a_sample_there_is_none():
    trex = TRex()
    assert trex.sample is None
    assert "sample" not in trex.components


def test_a_sample_goes_after_the_sample_position():
    trex = TRex(sample=aluminum())
    assert list(trex.components)[-5:] == [
        "sample_pos",
        "sample_monitor",
        "sample",
        "beamstop",
        "banana_detector",
    ]
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
    assert trex.arm_nboa.at == pytest.approx([-0.002022816656807591, 0, 1.895927], abs=1e-12)
    assert trex.arm_nboa.rotated == pytest.approx([0, -0.21499811064881996, 0], abs=1e-12)
    assert trex.arm_nboa.relative is None
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
        relative="arm_nboa",
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
    assert trex.guide_2.relative == trex.guide_3.relative == "arm_nboa"
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
    assert TRex().arm_nboa.rotated[1] == pytest.approx(-0.42999016682346453 / 2, abs=1e-5)


def test_the_center_line_tilts_toward_minus_x():
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
    instr_param_dict = {p.name: p.value for p in trex.instr_params}
    assert instr_param_dict["bender_in"] == 0
    assert instr_param_dict["b_rot"] == trex.params.b_rot


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
    components = [name for name, c in trex.components.items() if not isinstance(c, NXdisk_chopper)]
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
    names = [n for n in trex.components if n != "beamstop_monitor"]  # 3 m past the sample
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


def test_the_pre_sample_monitor_is_between_the_slit_and_the_sample():
    trex = TRex()
    assert trex.pre_sample_monitor.relative == "arm_mc"
    assert trex.pre_sample_monitor.at[2] == pytest.approx(163.5 - 161.985, abs=1e-12)
    assert trex.slit_mc.at[2] < trex.pre_sample_monitor.at[2] < trex.sample_pos.at[2]


def test_the_beamstop_is_just_in_front_of_the_banana_with_its_monitor_on_its_face():
    trex = TRex()
    monitor, stop, banana = trex.beamstop_monitor, trex.beamstop, trex.banana_detector
    positions = trex.positions()
    assert positions["beamstop_monitor"] == positions["beamstop"]
    assert (stop.xwidth, stop.yheight) == (monitor.xwidth, monitor.yheight)
    assert (monitor.xwidth, monitor.yheight) == (
        trex.pre_sample_monitor.xwidth,
        trex.pre_sample_monitor.yheight,
    )
    assert stop.relative == banana.relative == "sample_pos"
    assert banana.curved_radius - stop.at[2] == pytest.approx(0.001, abs=1e-12)


def test_the_beamstop_monitor_sees_the_beam_before_the_sample_and_gives_it_back_unchanged():
    trex = TRex(sample=aluminum())
    names = list(trex.components)
    assert names.index("beamstop_monitor") < names.index("sample_pos") < names.index("sample")
    assert trex.beamstop_monitor.restore_neutron == 1
    assert trex.beamstop_monitor.relative == "arm_mc"  # sample_pos comes after it


def test_the_banana_detector_is_around_the_sample_and_lists_its_pixels():
    trex = TRex(sample=aluminum())
    banana = trex.banana_detector
    positions = trex.positions()
    x, y, z = np.subtract(positions["banana_detector"], positions["sample"])
    assert (x, z) == pytest.approx((0, 0), abs=1e-9)  # its axis is through the sample
    # lowered: seen from the sample, 25 deg below the beam to 15 deg above
    bottom, top = y - banana.yheight / 2, y + banana.yheight / 2
    assert math.degrees(math.atan2(bottom, 3.0)) == pytest.approx(-25, abs=0.1)
    assert math.degrees(math.atan2(top, 3.0)) == pytest.approx(15, abs=0.1)
    assert banana.curved_radius == 3.0
    assert banana.bin_limits["theta"] == (-36.0, 144.0, 288)
    assert "mantid banana" in banana.options
    assert "list all neutrons" in banana.options


def test_two_instruments_share_no_component():
    first, second = TRex(), TRex()
    assert first == second
    assert first.guide_1 is not second.guide_1
    assert first.bw1 is not second.bw1
    first.guide_1.w1 = 0.1
    assert second.guide_1.w1 == pytest.approx(0.089955)


def test_running_conditions_give_the_settings_and_the_chopper_frequencies():
    trex = TRex()
    assert trex.running_conditions() == (
        "TRex · FF · λ0 = 2 Å · RRM = 18 · PS 189/-189 Hz · MC 252/252 Hz"
    )
    trex.update_params(mode="RR", lambda0=4.0, rrm=12)
    assert trex.running_conditions() == (
        "TRex · RR · λ0 = 4 Å · RRM = 12 · PS 126/-126 Hz · MC 168/168 Hz"
    )
