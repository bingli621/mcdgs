from pathlib import Path

import h5py
import pytest

from mcdgs.instruments.trex.model import TRex
from mcdgs.mcstas import simulation
from mcdgs.mcstas.simulation import Simulation, _run_points_in_parallel


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
    sim = Simulation(trex, output_dir=tmp_path)
    with pytest.raises(ValueError, match="sample_pos"):
        getattr(sim, method)(**kwargs)
    assert not any(tmp_path.iterdir())  # no run folder


def test_a_run_from_the_sample_without_a_beam_file_raises_before_anything(tmp_path):
    sim = Simulation(TRex(), output_dir=tmp_path)
    with pytest.raises(FileNotFoundError, match="run_to_sample"):
        sim.run_from_sample()
    assert not any(tmp_path.iterdir())  # no run folder


def test_the_beam_file_found_is_the_one_given_else_the_simulations(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    sim = Simulation(TRex(), output_dir=tmp_path)
    with pytest.raises(FileNotFoundError):
        sim.find_beam_file()  # none given, none found
    sim.beam_file = tmp_path / "a" / "beam.mcpl.gz"
    assert sim.find_beam_file() == sim.beam_file
    assert sim.find_beam_file("b/beam.mcpl.gz") == tmp_path.resolve() / "b" / "beam.mcpl.gz"


def test_the_folders_are_in_the_work_folder(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    trex = TRex()
    sim = Simulation(trex)
    assert sim.output_dir == tmp_path / "mcstas_output"
    assert sim.input_dir == tmp_path / "mcstas_input" / trex.name
    sim = Simulation(trex, work_dir=tmp_path / "w", output_dir="out")
    assert sim.output_dir == Path("out")
    assert sim.input_dir == tmp_path / "w" / "mcstas_input" / trex.name


def test_the_beam_file_is_the_newest_one_found(tmp_path, finished_run):
    assert Simulation(TRex(), output_dir=tmp_path).beam_file is None
    finished_run(tmp_path / "old_beam", 1, beam=True)
    new = finished_run(tmp_path / "new_beam", 2, beam=True)
    finished_run(tmp_path / "newer_scan", 3)
    assert Simulation(TRex(), output_dir=tmp_path).beam_file == new / "beam.mcpl.gz"


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
    beam = finished_run(tmp_path / "1_run_to_sample", 1, beam=True)
    scans = [finished_run(tmp_path / f"{n}_scan_b_rot", n) for n in (2, 3)]
    sim = Simulation(TRex(), output_dir=tmp_path)
    assert sim.list_runs() == ["3_scan_b_rot", "2_scan_b_rot", "1_run_to_sample"]
    assert sim.latest_run("*_run_to_sample") == beam
    assert sim.clean_runs(keep=1) == [scans[0]]
    assert sim.list_runs() == ["3_scan_b_rot", "1_run_to_sample"]
