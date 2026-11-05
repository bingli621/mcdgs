import h5py
import numpy as np
import pytest

from mcdgs.core.components import Monitor_nD
from mcdgs.diagnostics.trex.bender import (
    WAVELENGTH_EDGES,
    BenderScan,
    read_spectra,
    transmission,
)
from mcdgs.diagnostics.trex.extraction import Extraction
from mcdgs.mcstas.builder import build_mcstas


def test_the_bender_scan_adds_a_wavelength_monitor_after_guide_4b():
    scan = BenderScan()
    *up_to_guide_4, last = Extraction().components
    assert last == "beam_monitor"
    assert list(scan.components) == [*up_to_guide_4, "monitor_after_4b", "beam_monitor"]
    monitor = scan.monitor_after_4b
    assert isinstance(monitor, Monitor_nD)
    assert monitor.relative == "guide_4b"
    assert monitor.options == "lambda limits=[0.5 4.9] bins=11"
    assert len(WAVELENGTH_EDGES) - 1 == 11


def test_the_bender_scan_builds_in_mcstas(tmp_path):
    inst = build_mcstas(BenderScan(), input_path=tmp_path)
    assert [c.name for c in inst.component_list][-2:] == ["monitor_after_4b", "beam_monitor"]


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


def _write_scan(path, spectra):
    """A file as McStas writes a scan, with one entry per (data, errors) in spectra."""
    wavelengths = np.array([0.7, 1.1, 1.5])
    with h5py.File(path, "w") as file:
        for n, (data, errors) in enumerate(spectra, start=1):
            components = file.create_group(f"entry{n}/instrument/components")
            components.create_group("0009_guide_4b")
            spectrum = components.create_group(
                "0010_monitor_after_4b/output/monitor_after_guide_4b_L"
            )
            spectrum["data"] = data
            spectrum["errors"] = errors
            spectrum["Wavelength__Angs_"] = wavelengths
    return wavelengths


def test_the_spectra_are_read_in_the_order_of_the_scan(tmp_path):
    # 11 entries: by name, entry10 and entry11 would come before entry2
    spectra = [(np.full(3, float(n)), np.full(3, n / 10)) for n in range(1, 12)]
    wavelengths = _write_scan(tmp_path / "mccode.h5", spectra)
    centers, intensity, error = read_spectra(tmp_path / "mccode.h5")
    np.testing.assert_array_equal(centers, wavelengths)
    assert intensity.shape == error.shape == (11, 3)
    np.testing.assert_array_equal(intensity[:, 0], np.arange(1, 12))
    np.testing.assert_allclose(error[:, 0], np.arange(1, 12) / 10)
