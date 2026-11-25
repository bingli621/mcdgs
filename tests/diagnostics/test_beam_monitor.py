import numpy as np
import pytest
import scipp as sc

from mcdgs.diagnostics.beam_monitor import count_rate

BINS = {
    "toa": sc.array(dims=["toa"], values=[0.0, 1.0, 2.0], unit="ms"),
    "wavelength": sc.array(dims=["wavelength"], values=[1.0, 2.0, 4.0], unit="Å"),
}


def test_the_rates_of_the_flux():
    intensity = np.array([[1.0, 2.0], [3.0, 4.0]])  # n/s per bin
    rate = count_rate(intensity, BINS, frequency=10.0)
    # 10 n/s over the pulses; 7 n/s in the second time bin, over 10 pulses of 1 ms each
    assert rate.average == pytest.approx(10.0)
    assert rate.peak == pytest.approx(7.0 / 10.0 / 1e-3)
    assert rate.peak_time == pytest.approx(1.5)


def test_the_efficiency_is_taken_at_the_centers_of_the_wavelength_bins():
    intensity = np.array([[1.0, 0.0], [0.0, 1.0]])
    rate = count_rate(intensity, BINS, lambda wavelength: wavelength, frequency=10.0)
    assert rate.average == pytest.approx(1.5 + 3.0)
    assert rate.peak_time == pytest.approx(1.5)
