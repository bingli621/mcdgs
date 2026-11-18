import math

import pytest

from mcdgs.core import INSTRUMENTS, speed_from_wavelength
from mcdgs.core.components import ESS_butterfly, NXdisk_chopper, Polygon_ESS_butterfly
from mcdgs.instruments import CSPEC
from mcdgs.instruments.cspec import PARAMS_FILE
from mcdgs.instruments.cspec.choppers import DISKS, chopper_timings
from mcdgs.instruments.cspec.model import Params
from mcdgs.instruments.cspec.sources import wavelength_band
from mcdgs.io import load, read_yaml, to_dict
from mcdgs.samples import aluminum


def _close(a, b):
    """`a` and `b` are the same nested dict, numbers equal to 1e-12 relative."""
    if isinstance(a, dict):
        return isinstance(b, dict) and list(a) == list(b) and all(_close(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return isinstance(b, list) and len(a) == len(b) and all(map(_close, a, b))
    numbers = (int, float)
    if isinstance(a, numbers) and isinstance(b, numbers) and bool not in (type(a), type(b)):
        return math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-15)
    return a == b and type(a) is type(b)


def test_cspec_is_an_instrument_a_file_can_name():
    assert INSTRUMENTS["CSPEC"] is CSPEC


def test_params_yaml_is_the_default_cspec_saved():
    # if this fails after a change to the model: save(CSPEC(), PARAMS_FILE) and review the diff
    assert _close(read_yaml(PARAMS_FILE), to_dict(CSPEC()))
    assert _close(to_dict(load(PARAMS_FILE)), to_dict(CSPEC()))


def test_the_source_comes_first_then_the_frame_of_the_sheet():
    assert list(CSPEC().components)[:4] == ["origin", "source", "arm_iscs", "guide_01_011"]


def test_the_choppers_are_in_beam_order_each_after_the_arm_of_its_gap():
    cspec = CSPEC()
    names = list(cspec.components)
    assert list(cspec.choppers) == ["bw1", "bw2", "bw3", "ps1", "ps2", "rrm", "mc1", "mc2"]
    for name, disk in DISKS.items():
        assert cspec.components[name].relative == disk.gap_arm
        assert names.index(disk.gap_arm) < names.index(name)


def test_the_sample_position_comes_last_without_a_sample():
    cspec = CSPEC()
    assert list(cspec.components)[-1] == "sample_pos"
    assert cspec.sample is None


def test_a_sample_goes_after_the_sample_position():
    cspec = CSPEC(sample=aluminum())
    assert list(cspec.components)[-2:] == ["sample_pos", "sample"]
    assert cspec.sample.relative == "sample_pos"


def test_a_sample_relative_to_another_component_is_refused():
    sample = aluminum()
    sample.relative = "arm_mc"
    with pytest.raises(ValueError, match="arm_mc"):
        CSPEC(sample=sample)


def test_each_chopper_window_is_as_wide_as_the_exit_of_the_guide_in_front():
    cspec = CSPEC()
    for name, guide in cspec.upstream_guides.items():
        assert cspec.components[name].xwidth == guide.w2


def test_the_frequencies_are_those_of_drcspec():
    params = Params(nu_mc=168.0, n=3)
    cspec = CSPEC(params=params)
    assert cspec.chopper_frequencies() == {
        "bw1": 14,
        "bw2": 14,
        "bw3": 14,
        "ps1": 84,
        "ps2": -84,
        "rrm": 56,
        "mc1": 168,
        "mc2": -168,
    }


@pytest.mark.parametrize("name", DISKS)
def test_the_first_slit_is_on_the_beam_when_the_central_wavelength_arrives(name):
    lambda0, t0_offset = 4.0, 0.002
    nu, delay = chopper_timings(lambda0, 168.0, 3, 14, t0_offset)[name]
    arrival = t0_offset + DISKS[name].distance / speed_from_wavelength(lambda0)
    # the angle the disk has turned past its zero mark at the arrival, in the sense of nu
    turned = 360 * nu * (arrival - delay)
    assert turned % 360 == pytest.approx(DISKS[name].phase() % 360)


def test_disks_with_the_axle_above_the_beam_are_at_180_and_below_at_0():
    above = {n for n, d in DISKS.items() if d.beam_angle == 180}
    assert above == {"bw1", "bw2", "bw3", "ps1", "mc1"}


def test_the_pulse_shaping_disks_have_three_slits_120_deg_apart():
    for name in ("ps1", "ps2"):
        assert DISKS[name].slits == ((0, 24.23), (120, 24.23), (240, 24.23))


def test_the_rrm_chopper_is_halfway_from_the_start_of_its_gap_to_mc1():
    assert DISKS["rrm"].z == pytest.approx(DISKS["mc1"].z / 2)
    assert DISKS["rrm"].distance < DISKS["mc1"].distance < DISKS["mc2"].distance


def test_each_chopper_takes_its_timing_from_the_instrument_parameters():
    cspec = CSPEC()
    values = {p.name: p.value for p in cspec.instr_params}
    for name, (nu, delay) in chopper_timings(3.0, 14.0, 1, 14, 0.0).items():
        assert cspec.components[name].nu == f"nu_{name}"
        assert values[f"nu_{name}"] == nu
        assert values[f"delay_{name}"] == delay


def test_without_lambda0_the_choppers_have_no_timing():
    cspec = CSPEC(params=Params(lambda0=None))
    with pytest.raises(ValueError, match="lambda0"):
        _ = cspec.instr_params


def test_by_default_the_source_is_not_masked():
    assert type(CSPEC().source) is ESS_butterfly


def test_a_masked_source_is_masked_by_the_eight_choppers():
    source = CSPEC(source_masked=True).source
    assert isinstance(source, Polygon_ESS_butterfly)
    assert len(source.choppers) == 8


@pytest.mark.parametrize("masked", [False, True])
def test_a_new_lambda0_moves_the_band_and_keeps_the_source_masked_or_not(masked):
    cspec = CSPEC(source_masked=masked)
    kind = type(cspec.source)
    cspec.update_params(lambda0=5.0)
    assert (cspec.source.Lmin, cspec.source.Lmax) == wavelength_band(5.0)
    assert type(cspec.source) is kind


def test_the_running_conditions_name_the_frequencies():
    cspec = CSPEC(params=Params(nu_mc=168.0, n=3))
    assert cspec.running_conditions() == (
        "CSPEC · λ0 = 3 Å · PS 84/-84 Hz · RRM 56 Hz · MC 168/-168 Hz"
    )


def test_every_chopper_is_a_disk_chopper():
    assert all(isinstance(c, NXdisk_chopper) for c in CSPEC().choppers.values())
