from pathlib import Path

import matplotlib.pyplot as plt
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

    def plot(self):
        plt.figure()  # as McStasScript draws, without returning it


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


def test_plot_returns_the_figure_mcstasscript_draws(tmp_path, fake_diagnostics):
    diag = McStasScriptDiagnostics(TRex(), Workspace(tmp_path))
    diag.add_point(after="source")
    diag.add_view("l")
    diag.run()
    fig = diag.plot()
    assert fig is plt.gcf()
    plt.close(fig)


def test_the_folders_are_those_of_the_simulation(tmp_path):
    diag = McStasScriptDiagnostics(TRex(), Workspace(tmp_path))
    assert diag.workspace is diag.simulation.workspace
    assert diag.workspace.root == tmp_path


def test_run_analysis_is_a_run_of_its_own_with_the_intensity_zoomed(tmp_path, monkeypatch):
    def show_diagram(inst, analysis):
        """McStasScript's show_diagram without McStas: the intensity and ray count axes."""
        assert analysis and inst._run_settings["NeXus"] is False
        Path(inst._run_settings["output_path"]).mkdir()  # as McStas does
        _, (_, intensity, rays) = plt.subplots(1, 3)
        intensity.set(xlabel="Intensity [n/s]", xlim=(1, 1e12))
        rays.set(xlabel="Ray count", xlim=(1, 1e6))
        plt.show()  # held back

    monkeypatch.setattr("mcstasscript.interface.instr.McStas_instr.show_diagram", show_diagram)
    monkeypatch.setattr(plt, "show", lambda: pytest.fail("shown"), raising=True)
    diag = McStasScriptDiagnostics(TRex(), Workspace(tmp_path))
    run, fig = diag.run_analysis("a", ncount=1e4, intensity_limits=(1e10, 1e11))
    assert run.run_type is RunType.ANALYSIS and (run.folder / "run.json").is_file()
    _, intensity, rays = fig.axes
    assert intensity.get_xlim() == (1e10, 1e11)
    assert rays.get_xlim() == pytest.approx((1e4, 1e5))  # the same ratio to the intensity
    plt.close(fig)
