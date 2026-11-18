import numpy as np
import pytest

from mcdgs.diagnostics.cspec.source import define_instrument
from mcdgs.instruments import CSPEC
from mcdgs.instruments.cspec.sources import wavelength_band


def test_the_monitor_is_the_entrance_of_the_guide_system_of_cspec():
    positions = define_instrument().positions()
    np.testing.assert_allclose(positions["source_monitor"], CSPEC().positions()["guide_01_011"])


def test_the_monitor_is_as_large_as_the_entrance():
    monitor = define_instrument().source_monitor
    guide = CSPEC().components["guide_01_011"]
    assert (monitor.xwidth, monitor.yheight) == pytest.approx((guide.w1, guide.h1))


def test_with_lambda0_the_source_has_the_band_of_cspec():
    source = define_instrument(lambda0=3.0).source
    assert (source.Lmin, source.Lmax) == wavelength_band(3.0)
