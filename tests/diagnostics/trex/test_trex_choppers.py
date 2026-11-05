import pytest

from mcdgs.diagnostics.trex.choppers import (
    ArrivalSummary,
    count_groups,
    main,
    modes_table,
    summarize,
)
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


def test_the_modes_table_has_a_line_per_mode():
    summary = ArrivalSummary(
        neutrons=1925,
        wavelength_min=1.234,
        wavelength_max=2.678,
        sub_pulses=13,
        first_arrival=60.4,
        last_arrival=110.6,
    )
    table = modes_table({"FF": summary, "RR": summary}, neutrons=1_000_000)
    lines = table.splitlines()
    assert lines[0].startswith("| Mode | Neutrons at the sample, of 1,000,000 |")
    assert lines[2] == '| `"FF"` | 1925 | 1.23 to 2.68 | 13 | 60 to 111 |'
    assert len(lines) == 4  # the header, its rule and the two modes


def test_main_writes_the_three_files_of_the_docs(tmp_path):
    paths = main(out=tmp_path / "figures", neutrons=20_000)
    assert [p.name for p in paths] == ["trex_print.txt", "tof_result.txt", "tof_modes.txt"]
    assert all(p.parent == tmp_path / "figures" for p in paths)
    assert "sample_monitor" in paths[0].read_text()  # print(trex) lists every component
    assert len(paths[2].read_text().splitlines()) == 4
