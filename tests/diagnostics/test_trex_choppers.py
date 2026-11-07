import matplotlib

matplotlib.use("Agg")

from mcdgs.diagnostics.trex_00_overview import write_printout
from mcdgs.diagnostics.trex_05_choppers import main
from mcdgs.instruments.trex.model import TRex


def test_main_draws_every_disk_and_the_acceptance_diagram(tmp_path):
    paths = main(out=tmp_path / "figures")
    assert [p.name for p in paths] == [
        "chopper_bw1.svg",
        "chopper_bw2.svg",
        "chopper_ps1.svg",
        "chopper_ps2.svg",
        "chopper_mc1.svg",
        "chopper_mc2.svg",
        "trex_acceptance.png",
    ]
    assert all(p.parent == tmp_path / "figures" and p.stat().st_size > 0 for p in paths)
    assert all(p.read_text().startswith("<svg") for p in paths[:-1])


def test_the_printout_lists_every_component(tmp_path):
    path = write_printout(TRex(), tmp_path)
    assert path.name == "trex_print.txt"
    assert "sample_monitor" in path.read_text()
