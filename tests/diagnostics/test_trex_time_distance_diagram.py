import matplotlib

matplotlib.use("Agg")

from mcdgs.diagnostics.trex_06_time_distance_diagram import main


def test_main_draws_the_diagram_and_lists_its_parameters(tmp_path):
    paths = main(out=tmp_path / "figures", neutrons=10_000)
    assert [p.name for p in paths] == ["trex_time_distance.png", "trex_time_distance_params.txt"]
    assert all(p.parent == tmp_path / "figures" and p.stat().st_size > 0 for p in paths)
    table = paths[1].read_text()
    assert "| `neutrons`, per pulse | 10,000 |" in table
    assert "| `n_pulses` | 2 |" in table
    assert "| `ps2` | 108.050 m | \u2212189 Hz |" in table
