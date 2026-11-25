import matplotlib

matplotlib.use("Agg")

import pytest

from mcdgs.diagnostics.cspec.time_distance import define_instrument, main
from mcdgs.instruments import CSPEC
from mcdgs.tof import beamline


def test_the_sample_monitor_is_at_the_sample_position_after_the_choppers():
    cspec = define_instrument()
    assert list(cspec.components)[-2:] == ["sample_pos", "sample_monitor"]
    assert cspec.positions()["sample_monitor"] == pytest.approx(CSPEC().positions()["sample_pos"])


def test_tof_reaches_the_sample_monitor_after_the_last_chopper():
    _, after_source = beamline(define_instrument())
    distances = {name: distance for name, _, distance in after_source}
    assert distances["sample_monitor"] == max(distances.values()) > distances["mc2"]


def test_main_draws_the_diagram_and_lists_its_parameters(tmp_path):
    paths = main(out=tmp_path / "figures", neutrons=10_000)
    assert [p.name for p in paths] == ["cspec_time_distance.png", "cspec_time_distance_params.txt"]
    assert all(p.parent == tmp_path / "figures" and p.stat().st_size > 0 for p in paths)
    table = paths[1].read_text()
    assert "| `neutrons`, per pulse | 10,000 |" in table
    assert "| `nu_mc` | 224 Hz |" in table
    assert "| `t0_offset` | 1.7 ms |" in table
    assert "| `mc2` | 158.511 m | 224 Hz |" in table
