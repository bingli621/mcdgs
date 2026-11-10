"""Real McStas runs of a small instrument and of T-Rex: compiled and simulated, nothing faked."""

import json
import shutil
from dataclasses import dataclass

import h5py
import numpy as np
import pytest

from mcdgs.core.base import Instrument
from mcdgs.core.components import (
    Arm,
    ESS_butterfly,
    Monitor_nD,
    NXdisk_chopper,
    Polygon_ESS_butterfly,
    Progress_bar,
)
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


def _detected(instrument, tmp_path, seed=1, ncount=1e5):
    """The intensity on the detector of a run of ``instrument``."""
    h5_file = Simulation(instrument, work_dir=tmp_path).scan(
        "width", [0.1], ncount=ncount, seed=seed
    )
    return read_pixels(h5_file.parent, "detector", "source", "sample_pos")[0][1].sum()


def test_polygon_ess_butterfly_without_a_disk_is_ess_butterfly(tmp_path):
    plain = _Small()
    masked = _Small()
    masked.source = Polygon_ESS_butterfly(
        Lmin=2, Lmax=4, dist=2, focus_xw=0.05, focus_yh=0.05, save_polygons=0
    )
    # the same rays, from the same seed: nothing is outside the region, so nothing is drawn again
    assert _detected(masked, tmp_path / "masked").value == pytest.approx(
        _detected(plain, tmp_path / "plain").value, rel=1e-12
    )


def _chopped(source):
    """_Small with ``source`` and, 1 m from it, a disk turning at 14 Hz, open around 6 ms."""
    small = _Small()
    small.source = source
    # a 3 Å neutron leaving at 2 ms flies 1 m in 0.758 ms, so the slit is centered on it
    small.chopper = NXdisk_chopper(
        slit_edges=[-10, 10], radius=0.35, yheight=0.05, nu=14, delay=0.002758, at=[0, 0, 1]
    )
    return small


def test_polygon_ess_butterfly_gives_the_intensity_of_the_choppers_alone(tmp_path):
    # aperture: the 5 cm by 5 cm window of the disk, from chopper-lib's formula
    reach = np.sqrt(0.35**2 - 0.025**2)
    aperture = float(2 * np.degrees(np.arctan2(0.025, reach - 0.05)))
    disk = {"speed": 14, "delay": 0.002758, "beam": 0, "edges": [-10, 10], "path": 1.0, "aperture": aperture}  # fmt: skip
    shape = {"Lmin": 2, "Lmax": 4, "dist": 2, "focus_xw": 0.05, "focus_yh": 0.05}
    masked = _chopped(Polygon_ESS_butterfly(**shape, choppers=[disk]))

    h5_file = Simulation(masked, work_dir=tmp_path / "region").run(ncount=1e3, seed=1, mpi=None)
    region = json.loads((h5_file.parent / "source.json").read_text())  # named after the source
    # the disk passes a little over half of the 8.6 ms the source emits in, at each wavelength
    assert 0.5 < region["acceptance"] < 0.6

    # as much intensity from a quarter of the rays: over seeds, the two differ by 0.6 %,
    # each by about as much
    plain = _detected(_chopped(ESS_butterfly(**shape)), tmp_path / "plain", ncount=4e5)
    assert _detected(masked, tmp_path / "masked").value == pytest.approx(plain.value, rel=0.03)


@pytest.mark.parametrize("mode", ["FF", "RR"])
def test_trex_runs_through_with_its_choppers_turning(tmp_path, mode):
    # compiling alone let through an abort at INITIALIZE: the slit edges of the choppers
    trex = TRex()
    trex.update_params(mode=mode)
    h5_file = Simulation(trex, work_dir=tmp_path).run(ncount=1e3, seed=1, mpi=None)
    assert h5_file.is_file()
