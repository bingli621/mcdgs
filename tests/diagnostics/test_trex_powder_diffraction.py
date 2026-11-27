import pytest

from mcdgs.diagnostics.trex import diffraction
from mcdgs.diagnostics.trex.powder_diffraction import SAMPLE_SPLIT, define_instrument
from mcdgs.diagnostics.trex.solid_angles import detector_window


def test_the_instrument_has_the_powder_before_the_banana_detector():
    trex = define_instrument()
    names = list(trex.components)
    assert names.index("sample_pos") < names.index("sample") < names.index("banana_detector")
    assert type(trex.sample).__name__ == "PowderN"
    assert trex.sample.reflections == "Al2O3_sapphire.laz"
    assert trex.sample.split == SAMPLE_SPLIT
    assert type(trex.source).__name__ == "Polygon_ESS_butterfly"


def test_the_powder_spreads_its_cones_over_the_height_of_the_detector():
    assert define_instrument().sample.d_phi == pytest.approx(detector_window()["focus_ah"])


def test_the_slit_is_closed_to_the_size_of_the_powder():
    trex = define_instrument()
    assert trex.slit_mc.xwidth == pytest.approx(2 * trex.sample.radius)
    assert trex.slit_mc.yheight == pytest.approx(trex.sample.yheight)
    assert diffraction.define_instrument().slit_mc.xwidth > trex.slit_mc.xwidth  # the guide's
