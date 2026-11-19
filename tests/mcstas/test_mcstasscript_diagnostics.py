from pathlib import Path

import pytest

from mcdgs.instruments.trex.model import TRex
from mcdgs.mcstas import McStasScriptDiagnostics, RunType, Workspace


class FakeDiagnostics:
    """McStasScript's Diagnostics without McStas: it records what it is given."""

    def __init__(self, inst):
        self.inst, self.points, self.views, self.ran = inst, [], [], False

    def add_point(self, **point):
        self.points.append(point)

    def add_view(self, *axes, **kwargs):
        self.views.append((axes, kwargs))

    def run(self):
        Path(self.inst._run_settings["output_path"]).mkdir()  # as McStas does
        self.ran = True


@pytest.fixture
def fake_diagnostics(monkeypatch):
    monkeypatch.setattr("mcdgs.mcstas.mcstasscript_diagnostics.ms.Diagnostics", FakeDiagnostics)


def test_a_run_gives_the_points_and_views_to_mcstasscript_and_records_it(
    tmp_path, fake_diagnostics
):
    diag = McStasScriptDiagnostics(TRex(), Workspace(tmp_path))
    diag.add_point(after="source")
    diag.add_view("x", "y", bins=50)
    diag.add_view("l")
    run = diag.run("d", ncount=1e4, seed=3)
    fake = diag.diagnostics
    assert fake.ran and fake.points == [{"before": None, "after": "source", "rays": "all"}]
    assert fake.views == [(("x", "y"), {"colormap": "turbo", "bins": 50}), (("l", None), {})]
    assert run.run_type is RunType.DIAGNOSE and run.finished is not None
    assert run.folder == tmp_path / "runs" / "d" and (run.folder / "run.json").is_file()
    assert fake.inst._run_settings["NeXus"] is False  # McStasScript keeps them private
    assert (fake.inst._run_settings["ncount"], fake.inst._run_settings["seed"]) == (1e4, 3)


def test_every_run_gets_new_diagnostics_with_the_same_points_and_views(tmp_path, fake_diagnostics):
    diag = McStasScriptDiagnostics(TRex(), Workspace(tmp_path))
    diag.add_point(before="source")
    diag.add_view("l", colormap="jet")
    first = (diag.run("a"), diag.diagnostics)
    second = (diag.run("b", bender_in=1), diag.diagnostics)
    assert first[1] is not second[1] and first[1].views == second[1].views
    assert second[0].instr_params["bender_in"] == 1.0


def test_a_point_at_an_unknown_component_raises():
    diag = McStasScriptDiagnostics(TRex())
    with pytest.raises(KeyError, match="no_such_component"):
        diag.add_point(after="no_such_component")
    with pytest.raises(ValueError):
        diag.add_point()


def test_nothing_runs_or_plots_without_points_and_views(tmp_path):
    diag = McStasScriptDiagnostics(TRex(), Workspace(tmp_path))
    with pytest.raises(RuntimeError):
        diag.plot()
    diag.add_point(after="source")
    with pytest.raises(ValueError, match="view"):
        diag.run()
    assert not any(tmp_path.iterdir())  # no run folder
