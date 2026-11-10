"""Real McStas runs of a small instrument and of T-Rex: compiled and simulated, nothing faked."""

import shutil
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import pytest

from mcdgs.core.base import Instrument
from mcdgs.core.components import Arm, ESS_butterfly, Monitor_nD, Progress_bar
from mcdgs.core.params import InstrParam, param
from mcdgs.instruments.trex.model import TRex
from mcdgs.mcstas.nexus_io import read_pixels
from mcdgs.mcstas.simulation import Simulation

pytestmark = pytest.mark.skipif(shutil.which("mcrun") is None, reason="McStas is not on the PATH")

DETECTOR_AT = 1.0  # m from sample_pos


@dataclass(frozen=True)
class _Params:
    width: float = param(0.1, "m", "width of the detector")


class _Small(Instrument):
    params_class = _Params
    name = "small"

    def __init__(self, params=None):
        super().__init__(params)
        self.origin = Progress_bar()
        self.source = ESS_butterfly(Lmin=2, Lmax=4, dist=2, focus_xw=0.05, focus_yh=0.05)
        self.sample_pos = Arm(at=[0, 0, 2], rotated=[0, 0.5, 0])  # the detector stays in the beam
        self.detector = Monitor_nD(
            xwidth="width",
            yheight=0.1,
            options="x bins=4, y bins=3, list all, pixel id, t",
            at=[0, 0, DETECTOR_AT],
            relative="sample_pos",
        )

    @property
    def instr_params(self):
        return [InstrParam("width", self.params.width, "m", "width of the detector")]


@pytest.fixture
def sim(tmp_path):
    return Simulation(_Small(), work_dir=tmp_path)


def test_a_run_of_mcstas_writes_the_file_it_returns(sim):
    h5_file = sim.run(ncount=1e4, seed=1, mpi=None)
    assert h5_file.is_file()
    assert sim.list_runs() == [h5_file.parent.name]


def test_a_scan_runs_each_point_with_its_own_values_and_seed(sim):
    h5_file = sim.scan("width", [0.1, 0.05], ncount=1e4, seed=1, n_jobs=2)
    points = read_pixels(h5_file.parent, "detector", "source", "sample_pos")
    assert [params["width"] for params, _ in points] == [0.1, 0.05]
    with h5py.File(h5_file) as f:
        assert [f[f"entry{n}/simulation"].attrs["Seed"] for n in (1, 2)] == [b"2", b"3"]
    narrow, wide = points[1][1].sum().value, points[0][1].sum().value
    assert narrow < wide  # the instrument parameter reached McStas
    assert not (h5_file.parent / "point_0000").exists()  # merged, then removed


def test_a_split_run_gives_the_intensity_and_geometry_of_a_direct_run(sim):
    h5_file = sim.scan("width", [0.1], ncount=1e5, seed=1, n_jobs=1)
    direct = read_pixels(h5_file.parent, "detector", "source", "sample_pos")[0][1]

    beam = sim.run_to_sample(ncount=1e5, mpi=None)
    assert beam == sim.beam_file and beam.is_file()
    h5_file = sim.run_from_sample(seed=2, mpi=None)
    split = read_pixels(h5_file.parent, "detector", "source", "sample_pos", [0, 0, 1])[0][1]

    assert split.sum().value == pytest.approx(direct.sum().value, rel=0.05)
    # the part from the sample runs in the frame of sample_pos, at its origin
    assert np.allclose(split.coords["sample_position"].value, [0, 0, 0])
    pixel = split.coords["position"].values - split.coords["sample_position"].value
    assert np.allclose(pixel[:, 2], DETECTOR_AT)

    h5_file = sim.scan_from_sample("width", [0.1, 0.05], seed=3, n_jobs=2)
    points = read_pixels(h5_file.parent, "detector", "source", "sample_pos", [0, 0, 1])
    assert [params["width"] for params, _ in points] == [0.1, 0.05]
    assert points[1][1].sum().value < points[0][1].sum().value


@pytest.mark.parametrize("mode", ["FF", "RR"])
def test_trex_runs_through_with_its_choppers_turning(tmp_path, mode):
    # compiling alone let through an abort at INITIALIZE: the slit edges of the choppers
    trex = TRex()
    trex.update_params(mode=mode)
    h5_file = Simulation(trex, work_dir=tmp_path).run(ncount=1e3, seed=1, mpi=None)
    assert h5_file.is_file()


def test_train_mcpl_split_runs_end_to_end(tmp_path):
    component = (
        Path(__file__).parents[2]
        / "src/mcdgs/instruments/trex/tof_train_components/MCPL_output.comp"
    )
    if not component.is_file():
        pytest.skip("experimental MCPL_output component not yet available")
    sim = Simulation(TRex(), work_dir=tmp_path, tof_train=True)
    beam = sim.run_to_sample(ncount=10000, seed=1, mpi=None)
    assert beam.is_file()
    assert sim.run_from_sample(seed=2, mpi=None).is_file()
    assert sim.scan_from_sample("b_rot", [0], seed=3, n_jobs=1).is_file()
    assert sim.tof_train is True
