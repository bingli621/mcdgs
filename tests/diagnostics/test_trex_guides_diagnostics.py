import matplotlib

matplotlib.use("Agg")

from types import SimpleNamespace

import matplotlib.pyplot as plt

from mcdgs.diagnostics.trex import guides_diagnostics
from mcdgs.mcstas import McStasScriptDiagnostics


def test_main_saves_the_diagnostics_and_the_diagram(tmp_path, monkeypatch):
    runs = []

    def run(self, folder_name, **kwargs):
        runs.append((folder_name, kwargs, list(self._points)))
        return SimpleNamespace(folder=tmp_path / folder_name)  # as a record

    def run_analysis(self, folder_name, **kwargs):
        runs.append((folder_name, kwargs))
        return None, plt.figure()

    monkeypatch.setattr(McStasScriptDiagnostics, "run", run)
    monkeypatch.setattr(McStasScriptDiagnostics, "plot", lambda self: plt.figure())
    monkeypatch.setattr(McStasScriptDiagnostics, "run_analysis", run_analysis)
    paths = guides_diagnostics.main(out=tmp_path, ncount=1e3, analysis_ncount=1e4, seed=2)
    assert [p.name for p in paths] == ["trex_guides_diagnostics.png", "trex_intensity.png"]
    assert all(p.stat().st_size > 0 for p in paths)
    (name, kwargs, points), analysis = runs
    assert (name, kwargs) == ("trex_guides_diagnostics", {"ncount": 1e3, "seed": 2})
    assert [p["after"] or p["before"] for p in points] == [
        "source",
        "bw1",
        "guide_29",
        "guide_50a",
        "sample_pos",
    ]
    assert analysis == ("trex_guides_diagnostics_analysis", {"ncount": 1e4, "seed": 2})


def test_the_instrument_is_the_white_beam_with_every_chopper_parked():
    trex = guides_diagnostics.define_instrument()
    assert trex.params.lambda0 is None
    assert (trex.source.Lmin, trex.source.Lmax) == (0.1, 10.0)
    choppers = {p.name: p.value for p in trex.instr_params if p.name.startswith(("nu_", "delay_"))}
    assert len(choppers) == 12 and set(choppers.values()) == {0.0}
