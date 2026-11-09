import numpy as np
import pytest
import scipp as sc

from mcdgs.reduction.patches import Patches, find_patches, subpulse_flight

ORIGIN = sc.vector([0.0, 0.0, 0.0], unit="m")
END = sc.vector([0.0, 0.0, 3.0], unit="m")
SPEED_PER_ANGSTROM = 3956.034  # m/s times Å: the speed of a neutron of 1 Å


def _patches(time: list[float], intensity: list[float], wavelength: float = 2.0) -> Patches:
    ones = np.ones(len(time))
    return Patches(np.array(time), ones, wavelength * ones, ones, np.array(intensity))


def test_a_gap_longer_than_the_limit_starts_a_new_patch():
    bins = {
        "toa": sc.linspace("toa", 0.0, 10.0, 11, unit="ms"),
        "wavelength": sc.linspace("wavelength", 1.0, 2.0, 2, unit="Å"),
    }
    intensity = np.array([[1], [0], [1], [0], [0], [0], [1], [0], [0], [0]])
    patches = find_patches(intensity, bins, gap=1.0)  # 1 ms gap kept, 3 ms gap splits
    assert patches.intensity.tolist() == [2.0, 1.0]
    assert patches.time.tolist() == pytest.approx([1.5, 6.5])


def test_an_empty_map_has_no_patches():
    bins = {
        "toa": sc.linspace("toa", 0.0, 2.0, 3, unit="ms"),
        "wavelength": sc.linspace("wavelength", 1.0, 2.0, 2, unit="Å"),
    }
    assert find_patches(np.zeros((2, 1)), bins).time.size == 0


def test_the_weak_patches_are_not_sub_pulses():
    strong = _patches([1.0, 2.0, 3.0], [100.0, 1.0, 2.0]).strongest(0.02)
    assert strong.time.tolist() == [1.0, 3.0]


def test_each_sub_pulse_flies_from_its_patch_at_the_start_to_its_patch_at_the_end():
    wavelengths = [2.0, 1.8]
    flight_times = [1e3 * 3.0 * w / SPEED_PER_ANGSTROM for w in wavelengths]
    flight = subpulse_flight(
        _patches([10.0, 14.0], [100.0, 80.0]),
        ORIGIN,
        _patches([10.0 + flight_times[0], 14.0 + flight_times[1]], [90.0, 70.0]),
        END,
    )
    assert flight.dims == ("subpulse",)
    assert flight.values.tolist() == [100.0, 80.0]  # the intensity at the start
    assert flight.coords["time_at_start"].values.tolist() == [10.0, 14.0]
    assert flight.coords["wavelength_at_start"].values.tolist() == [2.0, 2.0]
    assert flight.coords["incident_wavelength"].to(unit="Å").values == pytest.approx(
        wavelengths, rel=1e-6
    )


def test_a_different_number_of_sub_pulses_raises():
    with pytest.raises(ValueError, match="2 sub-pulses at the start, 1 at the end"):
        subpulse_flight(_patches([10.0, 14.0], [1.0, 1.0]), ORIGIN, _patches([11.0], [1.0]), END)


def test_a_sub_pulse_that_reaches_the_end_first_raises():
    with pytest.raises(ValueError, match="reaches the end first"):
        subpulse_flight(_patches([10.0], [1.0]), ORIGIN, _patches([9.0], [1.0]), END)
