import matplotlib

matplotlib.use("Agg")

from mcdgs.diagnostics.cspec.time_distance import define_instrument, main
from mcdgs.instruments import CSPEC
from mcdgs.tof import build_tof
from mcdgs.tof.builder import DETECTOR


def test_the_instrument_is_cspec_at_its_defaults():
    assert define_instrument() == CSPEC()


def test_tof_reaches_the_detector_beyond_the_sample_monitor_after_the_last_chopper():
    model = build_tof(define_instrument(), neutrons=1000, seed=1)
    distances = {name: d.distance.value for name, d in model.detectors.items()}
    last_chopper = model.choppers["mc2"].distance.value
    assert distances[DETECTOR] == max(distances.values()) > distances["sample_monitor"]
    assert distances["sample_monitor"] > last_chopper


def test_main_draws_the_diagram_and_lists_its_parameters(tmp_path):
    paths = main(out=tmp_path / "figures", neutrons=10_000)
    assert [p.name for p in paths] == ["cspec_time_distance.png", "cspec_time_distance_params.txt"]
    assert all(p.parent == tmp_path / "figures" and p.stat().st_size > 0 for p in paths)
    table = paths[1].read_text()
    assert "| `neutrons`, per pulse | 10,000 |" in table
    assert "| `nu_mc` | 224 Hz |" in table
    assert "| `t0_offset` | 1.7 ms |" in table
    assert "| `mc2` | 158.511 m | 224 Hz |" in table
