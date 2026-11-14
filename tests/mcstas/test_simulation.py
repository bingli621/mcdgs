from pathlib import Path

import h5py
import pytest

from mcdgs.instruments.trex.model import TRex
from mcdgs.mcstas import simulation
from mcdgs.mcstas.runs import RunType
from mcdgs.mcstas.simulation import Simulation, _run_points_in_parallel
from mcdgs.mcstas.workspace import Workspace


def test_a_scan_has_a_point_per_value_of_its_variable():
    points = Simulation(TRex()).get_scan_points("b_rot", range(3), {"bender_in": 1})
    assert points == [
        {"bender_in": 1, "b_rot": 0.0},
        {"bender_in": 1, "b_rot": 1.0},
        {"bender_in": 1, "b_rot": 2.0},
    ]


@pytest.mark.parametrize(("variable", "instr_params"), [("b_rto", {}), ("b_rot", {"bender": 1})])
def test_a_scan_of_an_unknown_parameter_raises(variable, instr_params):
    with pytest.raises(KeyError):
        Simulation(TRex()).get_scan_points(variable, [1, 2], instr_params)


@pytest.mark.parametrize(
    ("method", "kwargs"),
    [
        ("run_to_sample", {}),
        ("run_from_sample", {"beam_file": "beam.mcpl"}),
        ("scan_from_sample", {"variable": "b_rot", "values": [1], "beam_file": "beam.mcpl"}),
    ],
)
def test_a_split_run_without_sample_pos_raises_before_anything(tmp_path, method, kwargs):
    trex = TRex()
    del trex.sample_pos
    sim = Simulation(trex, Workspace(tmp_path))
    with pytest.raises(ValueError, match="sample_pos"):
        getattr(sim, method)(**kwargs)
    assert not any(tmp_path.iterdir())  # no run folder


def test_a_run_from_the_sample_without_a_beam_file_raises_before_anything(tmp_path):
    sim = Simulation(TRex(), Workspace(tmp_path))
    with pytest.raises(FileNotFoundError, match="run_to_sample"):
        sim.run_from_sample()
    assert not any(tmp_path.iterdir())  # no run folder


def test_the_beam_file_found_is_the_one_given_else_the_simulations(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    sim = Simulation(TRex(), Workspace(tmp_path))
    with pytest.raises(FileNotFoundError):
        sim.find_beam_file()  # none given, none found
    sim.beam_file = tmp_path / "a" / "beam.mcpl.gz"
    assert sim.find_beam_file() == sim.beam_file
    assert sim.find_beam_file("b/beam.mcpl.gz") == tmp_path.resolve() / "b" / "beam.mcpl.gz"


def test_the_folders_are_those_of_the_workspace(tmp_path, monkeypatch):
    monkeypatch.delenv("MCDGS_WORKSPACE", raising=False)
    monkeypatch.chdir(tmp_path)
    sim = Simulation(TRex())
    assert sim.workspace == Workspace(tmp_path / "mcdgs_workspace")
    monkeypatch.chdir(tmp_path.parent)  # fixed when the simulation was created
    assert sim.workspace.root == tmp_path / "mcdgs_workspace"
    assert Simulation(TRex(), Workspace(tmp_path / "w")).workspace.root == tmp_path / "w"


def test_the_beam_file_is_the_newest_one_found(tmp_path, finished_run):
    assert Simulation(TRex(), Workspace(tmp_path)).beam_file is None
    finished_run(tmp_path / "runs" / "old_beam", 1, beam=True)
    new = finished_run(tmp_path / "runs" / "new_beam", 2, beam=True)
    finished_run(tmp_path / "runs" / "newer_scan", 3)
    assert Simulation(TRex(), Workspace(tmp_path)).beam_file == new / "beam.mcpl.gz"


def test_the_record_of_a_run_has_every_instrument_parameter_it_used(tmp_path):
    sim = Simulation(TRex(), Workspace(tmp_path))
    run, _ = sim._prepare("mine", RunType.RUN, {"bender_in": 1, "omega": 30}, ncount=1e4, seed=3)
    defaults = {p.name: p.value for p in TRex().instr_params}
    assert run.instr_params.keys() == defaults.keys()
    assert (run.instr_params["bender_in"], run.instr_params["omega"]) == (1.0, 30.0)  # as given
    assert run.instr_params["b_rot"] == defaults["b_rot"]  # not given: as the instrument has it
    assert (run.instrument, run.ncount, run.seed) == ("TRex", 1e4, 3)
    assert run.folder == tmp_path / "runs" / "mine" and not run.folder.exists()  # not saved yet


def test_the_record_of_a_scan_has_its_variable_apart_from_the_other_parameters(tmp_path):
    sim = Simulation(TRex(), Workspace(tmp_path))
    run, _ = sim._prepare(
        None, RunType.SCAN_FROM_SAMPLE, {"bender_in": 1}, variable="omega", values=[0.0, 1.0]
    )
    assert run.folder.name.endswith("_scan_omega_from_sample")
    assert (run.variable, run.values) == ("omega", [0.0, 1.0])
    assert "omega" not in run.instr_params and run.instr_params["bender_in"] == 1.0
    assert run.started.tzinfo is not None


def test_the_beam_run_of_a_record_is_the_run_its_beam_file_is_in(tmp_path):
    sim = Simulation(TRex(), Workspace(tmp_path))
    beam_file = sim.find_beam_file(tmp_path / "runs" / "to_sample" / "beam.mcpl.gz")
    run, _ = sim._prepare("a", RunType.RUN_FROM_SAMPLE, {}, beam_file=beam_file)
    assert run.beam_run == "to_sample"
    elsewhere = sim.find_beam_file(tmp_path / "elsewhere" / "beam.mcpl.gz")
    run, _ = sim._prepare("b", RunType.RUN_FROM_SAMPLE, {}, beam_file=elsewhere)
    assert run.beam_run is None  # not in a run of this workspace


class _FakeRun:
    """A McStasScript instrument whose `backengine` writes a one-entry `mccode.h5` holding its
    values, or raises if its b_rot is `fail_at`; it must be compiled first, and not again."""

    def __init__(self, fail_at=None):
        self.fail_at, self.values, self.options, self.compiled = fail_at, {}, {}, 0

    def settings(self, **options):
        self.options.update(options)

    def set_parameters(self, **values):
        self.values.update(values)

    def backengine(self):
        assert self.compiled == 1 and self.options["force_compile"] is False
        if self.values["b_rot"] == self.fail_at:
            raise ValueError("Simulation failed and no data was written to disk")
        folder = Path(self.options["output_path"])
        folder.mkdir()
        with h5py.File(folder / "mccode.h5", "w") as f:
            f["entry1/b_rot"] = self.values["b_rot"]
            f["entry1/seed"] = self.options["seed"]


@pytest.fixture
def fake_compile(monkeypatch):
    """`_compile` counts on the instrument, instead of running mcrun."""

    def compile_(inst):
        inst.compiled += 1

    monkeypatch.setattr(simulation, "_compile", compile_)


@pytest.mark.usefixtures("fake_compile")
def test_points_run_in_parallel_merge_in_order_with_their_own_seeds(tmp_path):
    points = [{"b_rot": float(w)} for w in range(5)]
    merged = _run_points_in_parallel(_FakeRun(), tmp_path / "scan", points, seed=10, n_jobs=2)
    assert merged == tmp_path / "scan" / "mccode.h5"
    with h5py.File(merged) as f:
        assert [f[f"entry{n}/b_rot"][()] for n in range(1, 6)] == list(range(5))
        assert [f[f"entry{n}/seed"][()] for n in range(1, 6)] == list(range(11, 16))
    assert [p.name for p in (tmp_path / "scan").iterdir()] == ["mccode.h5"]


@pytest.mark.usefixtures("fake_compile")
def test_a_failed_point_raises_naming_it_and_keeps_the_folders(tmp_path):
    points = [{"b_rot": float(w)} for w in range(3)]
    with pytest.raises(RuntimeError, match=r"point 1 \{'b_rot': 1.0\}"):
        _run_points_in_parallel(_FakeRun(fail_at=1.0), tmp_path / "scan", points, n_jobs=1)
    assert not (tmp_path / "scan" / "mccode.h5").exists()
    assert (tmp_path / "scan" / "point_0000" / "mccode.h5").exists()


@pytest.mark.usefixtures("fake_compile")
def test_no_points_raises_before_anything(tmp_path):
    with pytest.raises(ValueError, match="at least one point"):
        _run_points_in_parallel(_FakeRun(), tmp_path / "scan", [])
    assert not (tmp_path / "scan").exists()


def test_a_simulation_finds_and_cleans_its_own_runs(tmp_path, finished_run):
    beam = finished_run(tmp_path / "runs" / "1_run_to_sample", 1, beam=True)
    scans = [finished_run(tmp_path / "runs" / f"{n}_scan_b_rot", n) for n in (2, 3)]
    sim = Simulation(TRex(), Workspace(tmp_path))
    assert sim.list_runs() == ["3_scan_b_rot", "2_scan_b_rot", "1_run_to_sample"]
    assert sim.latest_run("*_run_to_sample") == beam
    assert sim.clean_runs(keep=1) == [scans[0]]
    assert sim.list_runs() == ["3_scan_b_rot", "1_run_to_sample"]
