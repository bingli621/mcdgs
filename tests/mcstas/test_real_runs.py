"""Real McStas runs of a small instrument and of T-Rex: compiled and simulated, nothing faked."""

import json
import shutil
from dataclasses import dataclass

import h5py
import mcstasscript as ms
import numpy as np
import pytest
import scipp as sc

from mcdgs.core import InstrParam, Instrument, param
from mcdgs.core.components import (
    Arm,
    ESS_butterfly,
    Monitor_nD,
    NXdisk_chopper,
    Polygon_ESS_butterfly,
    Progress_bar,
)
from mcdgs.instruments.trex.model import TRex
from mcdgs.mcstas import (
    Run,
    RunType,
    Simulation,
    Workspace,
    read_mcstas_maps,
    read_pixels,
)

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
        # a histogram, after sample_pos: mcrun counts a value of a scan only if one prints its
        # "Detector:" line, which the event list of the detector does not
        self.monitor = Monitor_nD(
            xwidth=0.2,
            yheight=0.2,
            options="t limits=[0 0.1] bins=10",
            restore_neutron=1,
            at=[0, 0, DETECTOR_AT / 2],
            relative="sample_pos",
        )
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
    return Simulation(_Small(), Workspace(tmp_path))


def test_a_run_of_mcstas_writes_the_file_it_returns(sim):
    run = sim.run(ncount=1e4, seed=1, mpi=None)
    assert run.data_file.is_file()
    assert Run.load(run.folder) == run  # the record of the run, in its folder
    assert sim.workspace.list_runs() == [run]
    assert (run.instrument, run.run_type, run.ncount, run.seed) == ("small", RunType.RUN, 1e4, 1)
    assert run.instr_params == {"width": 0.1}  # not given: as the instrument has it
    assert run.finished > run.started and run.duration == run.finished - run.started


def test_a_scan_runs_each_point_with_its_own_values_and_seed(sim):
    run = sim.scan("width", [0.1, 0.05], ncount=1e4, seed=1, n_jobs=2)
    points = read_pixels(run.folder, "detector", "source", "sample_pos")
    assert [params["width"] for params, _ in points] == [0.1, 0.05]
    with h5py.File(run.data_file) as f:
        # mcrun's seeds: seed + n * 1024
        assert [f[f"entry{n}/simulation"].attrs["Seed"] for n in (1, 2)] == [b"1", b"1025"]
    narrow, wide = points[1][1].sum().value, points[0][1].sum().value
    assert narrow < wide  # the instrument parameter reached McStas
    assert not (run.folder / "0").exists()  # mcrun merged the values, then removed them
    assert Run.load(run.folder) == run
    assert (run.run_type, run.variable, run.values) == (RunType.SCAN, "width", [0.1, 0.05])
    assert run.instr_params == {} and run.seed == 1  # the variable is in the values


def test_the_maps_of_a_scan_are_read_per_entry_as_mcstasscript_reads_the_first(sim):
    run = sim.scan("width", [0.1, 0.05, 0.02], ncount=1e4, seed=1, n_jobs=3)
    bins = {"toa": sc.linspace("toa", 0.0, 0.1, 11, unit="s")}  # of the histogram, "t ..."
    maps = read_mcstas_maps(run, {"monitor": bins})["monitor"]
    assert maps.dims == ("width", "toa")
    assert sc.identical(maps.coords["width"], sc.array(dims=["width"], values=run.values, unit="m"))
    # McStasScript reads the first entry alone: the same map
    first = ms.name_search("monitor", ms.load_data(str(run.folder)))
    np.testing.assert_array_equal(maps["width", 0].values, first.Intensity)
    single = sim.run(ncount=1e4, seed=1, mpi=None)  # as McStasScript loads it, without "width"
    assert read_mcstas_maps(single, {"monitor": bins})["monitor"].dims == ("toa",)
    with pytest.raises(ValueError, match=r"monitor: McStas has a map of \(10,\) bins"):
        read_mcstas_maps(run, {"monitor": {"toa": sc.linspace("toa", 0.0, 0.1, 6, unit="s")}})[
            "monitor"
        ]


def test_a_split_run_gives_the_intensity_and_geometry_of_a_direct_run(sim):
    direct_run = sim.scan("width", [0.1], ncount=1e5, seed=1, n_jobs=1)
    direct = read_pixels(direct_run.folder, "detector", "source", "sample_pos")[0][1]

    to_sample = sim.run_to_sample(ncount=1e5, mpi=None)
    assert to_sample.beam_file.is_file() and sim.beam_run == to_sample.name
    assert to_sample.run_type is RunType.RUN_TO_SAMPLE and Run.load(to_sample.folder) == to_sample
    from_sample = sim.run_from_sample(seed=2, mpi=None)
    assert Run.load(from_sample.folder) == from_sample
    assert from_sample.run_type is RunType.RUN_FROM_SAMPLE
    assert from_sample.beam_run == to_sample.name and from_sample.ncount is None
    split = read_pixels(from_sample.folder, "detector", "source", "sample_pos", [0, 0, 1])[0][1]

    assert split.sum().value == pytest.approx(direct.sum().value, rel=0.05)
    # the part from the sample runs in the frame of sample_pos, at its origin
    assert np.allclose(split.coords["sample_position"].value, [0, 0, 0])
    pixel = split.coords["position"].values - split.coords["sample_position"].value
    assert np.allclose(pixel[:, 2], DETECTOR_AT)

    # an instrument parameter given to the part from the sample reaches McStas
    narrow = sim.run_from_sample(seed=2, mpi=None, width=0.05)
    assert narrow.instr_params == {"width": 0.05}
    pixels = read_pixels(narrow.folder, "detector", "source", "sample_pos", [0, 0, 1])[0][1]
    assert pixels.sum().value < 0.75 * split.sum().value

    scan = sim.scan_from_sample("width", [0.1, 0.05], seed=3, n_jobs=2)
    points = read_pixels(scan.folder, "detector", "source", "sample_pos", [0, 0, 1])
    assert [params["width"] for params, _ in points] == [0.1, 0.05]
    assert points[1][1].sum().value < points[0][1].sum().value
    assert (scan.run_type, scan.values) == (RunType.SCAN_FROM_SAMPLE, [0.1, 0.05])
    assert scan.beam_run == to_sample.name


def _detected(instrument, tmp_path, seed=1, ncount=1e5):
    """The intensity on the detector of a run of ``instrument``."""
    run = Simulation(instrument, Workspace(tmp_path)).scan("width", [0.1], ncount=ncount, seed=seed)
    return read_pixels(run.folder, "detector", "source", "sample_pos")[0][1].sum()


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

    run = Simulation(masked, Workspace(tmp_path / "region")).run(ncount=1e3, seed=1, mpi=None)
    region = json.loads((run.folder / "source.json").read_text())  # named after the source
    # the disk passes a little over half of the 8.6 ms the source emits in, at each wavelength
    assert 0.5 < region["acceptance"] < 0.6

    # as much intensity from a quarter of the rays: over seeds, the two differ by 0.6 %,
    # each by about as much
    plain = _detected(_chopped(ESS_butterfly(**shape)), tmp_path / "plain", ncount=4e5)
    assert _detected(masked, tmp_path / "masked").value == pytest.approx(plain.value, rel=0.03)


def test_trex_with_the_masked_source_gives_the_intensity_on_the_sample(tmp_path):
    # the masked source gives the sample all the intensity of the plain one, from fewer rays:
    # at the defaults, 2.64e7 +- 5 % from 2e5 rays, and 2.72e7 +- 17 % from 2e6
    masked = TRex(source_masked=True)
    run = Simulation(masked, Workspace(tmp_path)).run(ncount=2e5, seed=1, mpi=None)
    region = json.loads((run.folder / "source.json").read_text())
    assert 0 < region["acceptance"] < 0.02  # 42 polygons, about 1 %
    with h5py.File(run.data_file) as f:
        (name,) = [k for k in f["entry1/data"] if k.startswith("sample_monitor")]
        intensity = f[f"entry1/data/{name}/data"][...].sum()
    assert intensity == pytest.approx(2.7e7, rel=0.2)


def test_a_split_run_compiles_with_the_masked_source(tmp_path):
    # the run from the sample leaves out the source, and the type of its chopper train with it
    small = _Small()
    disk = {"speed": 14, "delay": 0.002758, "beam": 0, "edges": [-10, 10], "path": 1.0, "aperture": 9}  # fmt: skip
    small.source = Polygon_ESS_butterfly(
        Lmin=2, Lmax=4, dist=2, focus_xw=0.05, focus_yh=0.05, choppers=[disk]
    )
    sim = Simulation(small, Workspace(tmp_path))
    assert sim.run_to_sample(ncount=1e4, mpi=None).beam_file.is_file()
    assert sim.run_from_sample(seed=2, mpi=None).data_file.is_file()


@pytest.mark.parametrize("mode", ["FF", "RR"])
def test_trex_runs_through_with_its_choppers_turning(tmp_path, mode):
    # compiling alone let through an abort at INITIALIZE: the slit edges of the choppers
    trex = TRex()
    trex.update_params(mode=mode)
    run = Simulation(trex, Workspace(tmp_path)).run(ncount=1e3, seed=1, mpi=None)
    assert run.data_file.is_file()
