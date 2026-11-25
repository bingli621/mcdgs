import numpy as np
import pytest
import scipp as sc

from mcdgs.core.components.source import L_MAX, L_MIN
from mcdgs.diagnostics.cspec.beam_monitor import define_instrument, intensity_near
from mcdgs.instruments import CSPEC


def test_the_monitors_and_the_guides_are_where_they_are_in_cspec():
    positions = define_instrument().positions()
    cspec_positions = CSPEC().positions()
    names = ["source_monitor", "guide_01_011", "guide_02_01", "arm_gap_bbg", "beam_monitor"]
    for name in names:
        np.testing.assert_allclose(positions[name], cspec_positions[name])


def test_the_beam_monitor_is_the_last_component():
    assert list(define_instrument().components)[-1] == "beam_monitor"


def test_the_source_and_the_beam_monitor_have_the_widest_band():
    instrument = define_instrument()
    assert (instrument.source.Lmin, instrument.source.Lmax) == (L_MIN, L_MAX)
    assert instrument.beam_monitor.bin_limits["lambda"][:2] == (L_MIN, L_MAX)


def _map() -> sc.DataArray:
    """Two time bins and four wavelength bins, 0.1 Å wide, from 2.8 Å."""
    return sc.DataArray(
        sc.array(dims=["toa", "wavelength"], values=np.arange(8.0).reshape(2, 4), unit="1/s"),
        coords={"wavelength": sc.linspace("wavelength", 2.8, 3.2, 5, unit="Å")},
    )


def test_intensity_near_sums_the_bins_whose_centers_are_near():
    # centers 2.85, 2.95, 3.05, 3.15: within 0.05 Å of 3 Å are the second and the third
    assert intensity_near(_map(), 3.0) == pytest.approx(1 + 2 + 5 + 6)


def test_intensity_near_raises_outside_the_map():
    with pytest.raises(ValueError, match="no wavelength bin"):
        intensity_near(_map(), 5.0)
