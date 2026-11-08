from dataclasses import dataclass, replace

import numpy as np
import pytest
import scipp as sc
import tof

from mcdgs.core.base import Instrument
from mcdgs.core.components import (
    Arm,
    ESS_butterfly,
    Guide_gravity,
    Monitor_nD,
    Progress_bar,
)
from mcdgs.core.params import InstrParam
from mcdgs.core.utils import evaluate, speed_from_wavelength
from mcdgs.instruments.trex.choppers import DISK_DISTANCE
from mcdgs.instruments.trex.model import TRex
from mcdgs.samples import aluminum
from mcdgs.tof.builder import build_tof


@dataclass(frozen=True)
class _NoParams:
    pass


class _Beamline(Instrument):
    """A source, a bent guide, two monitors and a sample position: no choppers."""

    params_class = _NoParams
    name = "beamline"

    def __init__(self, params=None):
        super().__init__(params)
        self.origin = Progress_bar(at=[0, 0, 0])
        self.source = ESS_butterfly(Lmin="lmin", Lmax=6.0, n_pulses=2, at=[0, 0, 0])
        self.guide = Guide_gravity(w1=0.05, h1=0.05, l=10, at=[0, 0, 2], relative="source")
        self.axis = Arm(at=[0, 0, 12], rotated=[0, 10, 0], relative="source")
        self.monitor = Monitor_nD(xwidth=0.1, yheight=0.1, at=[0, 0, 3], relative="axis")
        self.sample_monitor = Monitor_nD(xwidth=0.1, yheight=0.1, at=[0, 0, 4], relative="axis")
        self.sample_pos = Arm(at=[0, 0, 4], relative="axis")

    @property
    def instr_params(self):
        return [InstrParam("lmin", 1.0, "Å", "shortest wavelength")]


def test_source_takes_the_band_and_pulses_of_the_moderator():
    model = build_tof(_Beamline(), neutrons=1000, seed=1)
    source = model.source
    assert source.pulses == 2
    assert source.neutrons == 1000
    wavelength = source.data.coords["wavelength"]
    assert wavelength.min().value >= 1.0  # Lmin from the instrument parameter
    assert wavelength.max().value <= 6.0


def test_source_stays_where_tof_sampled_its_pulse():
    source = build_tof(_Beamline(), neutrons=1000, seed=1).source
    assert sc.identical(source.distance, sc.scalar(0.05, unit="m"))  # 5 cm from the moderator


@pytest.mark.parametrize(
    ("sector", "beamline", "facility"), [("W", 7, "ess-trex"), ("N", 1, "ess")]
)
def test_source_is_the_pulse_of_the_beam_port(sector, beamline, facility):
    instrument = _Beamline()
    instrument.source = replace(instrument.source, sector=sector, beamline=beamline)
    assert build_tof(instrument, neutrons=1000, seed=1).source.facility == facility


def test_trex_gets_its_own_pulse():
    assert build_tof(TRex(), neutrons=1000, seed=1).source.facility == "ess-trex"


def test_monitors_become_detectors_and_the_rest_is_left_out():
    model = build_tof(_Beamline(), neutrons=1000, seed=1)
    assert list(model.components) == ["monitor", "sample_monitor"]
    for name, after_arm in (("monitor", 3), ("sample_monitor", 4)):
        # the arm, 12 m on, turns the beam by 10 deg: the straight line from the source
        x, z = after_arm * np.sin(np.radians(10)), 12 + after_arm * np.cos(np.radians(10))
        assert model.detectors[name].distance.value == pytest.approx(np.hypot(x, z)), name


def test_a_banana_detector_is_its_radius_beyond_the_sample():
    trex = TRex()
    model = build_tof(trex, neutrons=1000, seed=1)
    sample = model.detectors["sample_monitor"].distance.value
    banana = model.detectors["banana_detector"].distance.value
    # within 1 mm: the banana's center is 0.3 m below the sample, off the straight line
    assert banana - sample == pytest.approx(trex.banana_detector.curved_radius, abs=1e-3)
    assert banana > model.detectors["beamstop_detector"].distance.value


def test_the_model_runs():
    result = build_tof(_Beamline(), neutrons=1000, seed=1).run()
    assert isinstance(result, tof.Result)


def test_the_same_seed_gives_the_same_source():
    first, second = (build_tof(_Beamline(), neutrons=1000, seed=7) for _ in range(2))
    assert sc.identical(first.source.data, second.source.data)


def test_an_instrument_without_monitors_gives_no_detectors():
    beamline = _Beamline()
    beamline.monitor = beamline.sample_monitor = None
    assert not build_tof(beamline, neutrons=1000).detectors


@pytest.mark.parametrize("n_sources", [0, 2])
def test_the_instrument_needs_one_source(n_sources):
    beamline = _Beamline()
    if n_sources == 0:
        beamline.source = None
    else:
        beamline.second_source = ESS_butterfly(Lmin=1, Lmax=2)
    with pytest.raises(ValueError, match="tof needs one ESS_butterfly"):
        build_tof(beamline)


def test_a_component_tof_does_not_understand_is_left_out():
    beamline = _Beamline()
    beamline.sample = replace(aluminum(), relative="sample_pos")
    model = build_tof(beamline, neutrons=1000)
    assert list(model.detectors) == ["monitor", "sample_monitor"]
    assert not model.choppers


def _mcstas_open(chopper, params, t):
    """Whether the McStas disk lets a ray through on the beam axis at times `t` (s): the rule of
    NXdisk_chopper.comp's TRACE, written out again here, the reference for tof."""
    mark = 360 * evaluate(chopper.nu, params) * (t - evaluate(chopper.delay, params))
    on_disc = np.mod(evaluate(chopper.beam_angle, params) - mark, 360)
    edges = np.array(chopper.slit_edges, float)
    through = np.zeros(np.shape(t), bool)
    for opening, closing in zip(edges[0::2], edges[1::2], strict=True):
        through |= np.mod(on_disc - np.mod(opening, 360), 360) < closing - opening
    return through


def _tof_open(chopper, t):
    """Whether the tof chopper is open at times `t` (s)."""
    opens, closes = chopper.open_close_times(time_limit=sc.scalar(t.max(), unit="s"), unit="s")
    through = np.zeros(np.shape(t), bool)
    for opening, closing in zip(opens.values, closes.values, strict=True):
        through |= (t >= opening) & (t < closing)
    return through


def _trex(mode):
    trex = TRex()
    trex.update_params(mode=mode)
    return trex


def test_trex_choppers_are_placed_in_beam_order():
    model = build_tof(TRex(), neutrons=1000, seed=1)
    assert list(model.choppers) == list(DISK_DISTANCE)
    for name, chopper in model.choppers.items():
        assert chopper.distance.value == pytest.approx(DISK_DISTANCE[name], abs=1e-3), name
    assert list(model.detectors) == [
        "beam_monitor",
        "pre_sample_monitor",
        "sample_monitor",
        "beamstop_detector",
        "banana_detector",
    ]


@pytest.mark.parametrize("mode", ["FF", "RR"])
def test_tof_choppers_open_when_the_mcstas_disks_do(mode):
    # every 0.5 us over two pulses: under 2 % of the shortest opening, 28 us, that of the
    # monochromatic choppers
    trex = _trex(mode)
    params = trex.numeric_instr_params()
    t = np.arange(0, 2 / 14, 0.5e-6)
    for name, chopper in build_tof(trex, neutrons=1000, seed=1).choppers.items():
        mcstas = _mcstas_open(trex.components[name], params, t)
        assert mcstas.any(), name
        np.testing.assert_array_equal(_tof_open(chopper, t), mcstas, err_msg=name)


@pytest.mark.parametrize("mode", ["FF", "RR"])
def test_the_neutron_the_choppers_are_phased_to_gets_through_every_disk(mode):
    trex = _trex(mode)
    speed = speed_from_wavelength(trex.params.lambda0)
    for name, chopper in build_tof(trex, neutrons=1000, seed=1).choppers.items():
        arrival = trex.params.t0_offset + chopper.distance.value / speed
        assert _tof_open(chopper, np.array([arrival]))[0], name


def test_neutrons_reach_the_sample_exactly_when_every_mcstas_disk_lets_them_through():
    trex = TRex()
    params = trex.numeric_instr_params()
    model = build_tof(trex, neutrons=200_000, seed=1)
    neutrons = model.run()["sample_monitor"].data["pulse", 0]
    birth = neutrons.coords["birth_time"].to(unit="s").values
    speed = neutrons.coords["speed"].values
    start = model.source.distance.value  # the neutrons start where tof sampled the pulse
    through = np.ones(birth.shape, bool)
    for name, chopper in model.choppers.items():
        flight = (chopper.distance.value - start) / speed
        through &= _mcstas_open(trex.components[name], params, birth + flight)
    reached = ~neutrons.masks["blocked_by_others"].values
    assert reached.sum() > 100  # enough to mean something
    np.testing.assert_array_equal(reached, through)


def test_a_parked_disk_raises():
    trex = TRex()
    trex.bw1 = replace(trex.bw1, nu=0)
    with pytest.raises(ValueError, match="bw1: a parked disk"):
        build_tof(trex, neutrons=1000)


def test_trex_disks_lie_within_1_mm_of_the_distances_the_choppers_are_phased_to():
    trex = TRex()
    positions = trex.positions()
    source = np.array(positions["source"])
    for name, distance in DISK_DISTANCE.items():
        straight = np.linalg.norm(np.array(positions[name]) - source)
        assert straight == pytest.approx(distance, abs=1e-3), name
