import pytest

from mcdgs.diagnostics.trex.choppers import count_groups, summarize
from mcdgs.instruments.trex.model import TRex
from mcdgs.tof.builder import build_tof


def test_a_gap_longer_than_the_limit_starts_a_new_group():
    assert count_groups([5.0, 0.0, 10.0, 0.1, 5.2], gap=1.0) == 3  # order does not matter
    assert count_groups([0.0, 1.0, 2.0], gap=1.0) == 1  # a gap equal to the limit does not
    assert count_groups([3.0], gap=1.0) == 1


def test_there_must_be_times_to_group():
    with pytest.raises(ValueError, match="no times"):
        count_groups([], gap=1.0)


def test_trex_cuts_the_band_into_sub_pulses():
    result = build_tof(TRex(), neutrons=100_000, seed=1).run()
    summary = summarize(result, "sample_monitor")
    assert summary.neutrons > 0
    # the band the bandwidth choppers pass around lambda0 = 2 Å
    assert 1.3 < summary.wavelength_min < summary.wavelength_max < 2.7
    assert summary.sub_pulses == 13
    assert 55 < summary.first_arrival < summary.last_arrival < 115  # ms


def test_a_detector_no_neutron_reaches_raises():
    # about 1 neutron in 540 gets through: none of these 10
    result = build_tof(TRex(), neutrons=10, seed=1).run()
    with pytest.raises(ValueError, match="sample_monitor"):
        summarize(result, "sample_monitor")
