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
from mcdgs.instruments.trex.choppers import DISK_DISTANCE
from mcdgs.instruments.trex.model import TRex
from mcdgs.samples import aluminum
from mcdgs.tof.builder import build_tof


@dataclass(frozen=True)
class _NoParams:
    pass


class _Beamline(Instrument):
    """A source, a bent guide, a monitor and a sample position: no choppers."""

    params_class = _NoParams
    name = "beamline"

    def __init__(self, params=None):
        super().__init__(params)
        self.origin = Progress_bar(at=[0, 0, 0])
        self.source = ESS_butterfly(Lmin="lmin", Lmax=6.0, n_pulses=2, at=[0, 0, 0])
        self.guide = Guide_gravity(w1=0.05, h1=0.05, l=10, at=[0, 0, 2], relative="source")
        self.axis = Arm(at=[0, 0, 12], rotated=[0, 10, 0], relative="source")
        self.monitor = Monitor_nD(xwidth=0.1, yheight=0.1, at=[0, 0, 3], relative="axis")
        self.sample_pos = Arm(at=[0, 0, 4], relative="axis")

    @property
    def instr_params(self):
        return [InstrParam("lmin", 1.0, "Å", "shortest wavelength")]


def test_source_takes_the_band_and_pulses_of_the_moderator():
    model = build_tof(_Beamline(), neutrons=1000, seed=1)
    source = model.source
    assert source.pulses == 2
    assert source.neutrons == 1000
    assert sc.identical(source.distance, sc.scalar(0.0, unit="m"))
    wavelength = source.data.coords["wavelength"]
    assert wavelength.min().value >= 1.0  # Lmin from the instrument parameter
    assert wavelength.max().value <= 6.0


def test_monitors_and_named_components_become_detectors_and_optics_are_skipped():
    model = build_tof(_Beamline(), neutrons=1000, detectors=("sample_pos",), seed=1)
    assert list(model.components) == ["monitor", "sample_pos"]
    for name, after_arm in (("monitor", 3), ("sample_pos", 4)):
        # the arm, 12 m on, turns the beam by 10 deg: the straight line from the source
        x, z = after_arm * np.sin(np.radians(10)), 12 + after_arm * np.cos(np.radians(10))
        assert model.detectors[name].distance.value == pytest.approx(np.hypot(x, z)), name


def test_the_model_runs():
    result = build_tof(_Beamline(), neutrons=1000, detectors=("sample_pos",), seed=1).run()
    assert isinstance(result, tof.Result)


def test_the_same_seed_gives_the_same_source():
    first, second = (build_tof(_Beamline(), neutrons=1000, seed=7) for _ in range(2))
    assert sc.identical(first.source.data, second.source.data)


def test_a_detector_must_be_a_component_after_the_source():
    with pytest.raises(ValueError, match=r"no component \['origin', 'sample'\] after"):
        build_tof(_Beamline(), detectors=("origin", "sample"))


@pytest.mark.parametrize("n_sources", [0, 2])
def test_the_instrument_needs_one_source(n_sources):
    beamline = _Beamline()
    if n_sources == 0:
        beamline.source = None
    else:
        beamline.second_source = ESS_butterfly(Lmin=1, Lmax=2)
    with pytest.raises(ValueError, match="tof needs one ESS_butterfly"):
        build_tof(beamline)


def test_a_component_without_a_conversion_raises():
    beamline = _Beamline()
    beamline.sample = replace(aluminum(), relative="sample_pos")
    with pytest.raises(ValueError, match=r"sample \(Single_crystal\): no tof conversion"):
        build_tof(beamline, neutrons=1000)


def test_choppers_are_not_converted_yet():
    with pytest.raises(NotImplementedError, match="bw1: NXdisk_chopper"):
        build_tof(TRex(), neutrons=1000)


def test_trex_disks_lie_within_1_mm_of_the_distances_the_choppers_are_phased_to():
    trex = TRex()
    positions = trex.positions()
    source = np.array(positions["source"])
    for name, distance in DISK_DISTANCE.items():
        straight = np.linalg.norm(np.array(positions[name]) - source)
        assert straight == pytest.approx(distance, abs=1e-3), name
