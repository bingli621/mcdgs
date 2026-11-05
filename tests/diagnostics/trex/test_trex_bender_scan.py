import h5py
import numpy as np
import pytest

from mcdgs.diagnostics.trex.bender_scan import BAND_BINS, read_spectra, transmission


def test_the_transmission_divides_each_band_by_the_reference():
    intensity = np.array([[2.0, 3.0], [4.0, 6.0]])  # 2 scan points, 2 bands
    reference = np.array([4.0, 6.0])
    ratio, _ = transmission(intensity, np.zeros((2, 2)), reference, np.zeros(2))
    np.testing.assert_allclose(ratio, [[0.5, 0.5], [1.0, 1.0]])


def test_the_relative_errors_of_the_transmission_add_in_quadrature():
    # relative errors 3 % and 4 %: 5 % of the ratio
    ratio, ratio_error = transmission(
        np.array([[10.0]]), np.array([[0.3]]), np.array([20.0]), np.array([0.8])
    )
    assert ratio[0, 0] == pytest.approx(0.5)
    assert ratio_error[0, 0] == pytest.approx(0.05 * 0.5)


def _write_scan(path, maps):
    """A file as McStas writes a scan, with one entry per (data, errors) map of the beam monitor."""
    wavelengths = np.linspace(0.05, 1.55, 2 * BAND_BINS)  # bin centers: two bands
    with h5py.File(path, "w") as file:
        for n, (data, errors) in enumerate(maps, start=1):
            components = file.create_group(f"entry{n}/instrument/components")
            components.create_group("0009_guide_4b")
            time_map = components.create_group("0010_beam_monitor/output/beam_monitor_t_L")
            time_map["data"] = data
            time_map["errors"] = errors
            time_map["Wavelength__Angs_"] = wavelengths
    return wavelengths


def test_the_spectra_are_read_in_the_order_of_the_scan(tmp_path):
    # 11 entries: by name, entry10 and entry11 would come before entry2
    shape = (3, 2 * BAND_BINS)  # 3 time bins
    maps = [(np.full(shape, float(n)), np.full(shape, 0.0)) for n in range(1, 12)]
    _write_scan(tmp_path / "mccode.h5", maps)
    _, intensity, _ = read_spectra(tmp_path / "mccode.h5")
    assert intensity.shape == (11, 2)
    # each band sums 3 time bins and BAND_BINS wavelength bins
    np.testing.assert_allclose(intensity[:, 0], np.arange(1, 12) * 3 * BAND_BINS)


def test_a_band_sums_its_bins_and_their_errors_in_quadrature(tmp_path):
    shape = (4, 2 * BAND_BINS)  # 4 time bins
    wavelengths = _write_scan(tmp_path / "mccode.h5", [(np.ones(shape), np.full(shape, 0.5))])
    centers, intensity, error = read_spectra(tmp_path / "mccode.h5")
    np.testing.assert_allclose(centers, wavelengths.reshape(2, BAND_BINS).mean(axis=1))
    n_bins = 4 * BAND_BINS
    np.testing.assert_allclose(intensity, [[n_bins, n_bins]])
    np.testing.assert_allclose(error, [[0.5 * np.sqrt(n_bins)] * 2])


def test_bins_that_are_not_whole_bands_are_refused(tmp_path):
    shape = (1, BAND_BINS + 1)
    with h5py.File(tmp_path / "mccode.h5", "w") as file:
        time_map = file.create_group("entry1/instrument/components/0010_beam_monitor").create_group(
            "output/beam_monitor_t_L"
        )
        time_map["data"] = np.ones(shape)
        time_map["errors"] = np.ones(shape)
        time_map["Wavelength__Angs_"] = np.arange(BAND_BINS + 1)
    with pytest.raises(ValueError, match="whole number of bands"):
        read_spectra(tmp_path / "mccode.h5")
