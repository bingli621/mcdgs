from pathlib import Path

import h5py
import pytest

from mcdgs.mcstas import random_seed, run_points_in_parallel, scan


def test_a_random_seed_is_one_mcrun_accepts():
    seeds = [random_seed() for _ in range(20)]
    assert all(isinstance(seed, int) and 0 <= seed < 10**9 for seed in seeds)
    assert len(set(seeds)) > 1


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
    """`compile_instrument` counts on the instrument, instead of running mcrun."""

    def compile_(inst):
        inst.compiled += 1

    monkeypatch.setattr(scan, "compile_instrument", compile_)


@pytest.mark.usefixtures("fake_compile")
def test_points_run_in_parallel_merge_in_order_with_their_own_seeds(tmp_path):
    points = [{"b_rot": float(w)} for w in range(5)]
    merged = run_points_in_parallel(_FakeRun(), tmp_path / "scan", points, seed=10, n_jobs=2)
    assert merged == tmp_path / "scan" / "mccode.h5"
    with h5py.File(merged) as f:
        assert [f[f"entry{n}/b_rot"][()] for n in range(1, 6)] == list(range(5))
        assert [f[f"entry{n}/seed"][()] for n in range(1, 6)] == list(range(11, 16))
    assert [p.name for p in (tmp_path / "scan").iterdir()] == ["mccode.h5"]


@pytest.mark.usefixtures("fake_compile")
def test_a_failed_point_raises_naming_it_and_keeps_the_folders(tmp_path):
    points = [{"b_rot": float(w)} for w in range(3)]
    with pytest.raises(RuntimeError, match=r"point 1 \{'b_rot': 1.0\}"):
        run_points_in_parallel(_FakeRun(fail_at=1.0), tmp_path / "scan", points, n_jobs=1)
    assert not (tmp_path / "scan" / "mccode.h5").exists()
    assert (tmp_path / "scan" / "point_0000" / "mccode.h5").exists()


@pytest.mark.usefixtures("fake_compile")
def test_no_points_raises_before_anything(tmp_path):
    with pytest.raises(ValueError, match="at least one point"):
        run_points_in_parallel(_FakeRun(), tmp_path / "scan", [])
    assert not (tmp_path / "scan").exists()
