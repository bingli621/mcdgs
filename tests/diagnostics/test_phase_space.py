import numpy as np
import pytest

from mcdgs.diagnostics import phase_space as ps


def uniform_rays(n, x=(-0.07, 0.07), divergence=2.0, seed=0):
    """Rays of weight 1 spread evenly over x, |y| < 5 mm, ±divergence and 0.1 to 10 Å."""
    rng = np.random.default_rng(seed)
    return {
        "x": rng.uniform(*x, n),
        "y": rng.uniform(-0.005, 0.005, n),
        "dx": rng.uniform(-divergence, divergence, n),
        "dy": rng.uniform(-divergence, divergence, n),
        "l": rng.uniform(0.1, 10.0, n),
        "p": np.ones(n),
    }


def test_a_window_holds_the_rays_in_its_area_and_its_divergence():
    window = ps.Window(half_width=0.01, half_height=0.005, half_divergence=0.5, x_center=0.05)
    rays = {
        "x": np.array([0.05, 0.065, 0.05, 0.05]),
        "y": np.zeros(4),
        "dx": np.array([0.0, 0.0, 0.6, 0.0]),
        "dy": np.array([0.0, 0.0, 0.0, -0.6]),
    }
    assert window.contains(rays).tolist() == [True, False, False, False]
    assert window.area == pytest.approx(2.0)  # cm²: 2 × 1 cm


def test_the_intensity_in_a_window_is_the_weight_of_its_rays_in_each_bin():
    rays = uniform_rays(10_000)
    rays["p"] = np.full(10_000, 0.5)
    window = ps.Window(0.01, 0.005, 1.0)
    intensity, counts = ps.window_intensity(rays, window, np.array([0.1, 5.0, 10.0]))
    inside = window.contains(rays)
    assert counts.sum() == inside.sum()
    assert intensity.tolist() == (0.5 * counts).tolist()


def test_the_brilliance_of_an_even_beam_is_its_intensity_over_its_phase_space():
    # 2 cm × 1 cm, ±2° both ways, 9.9 Å: every window inside it has the same brilliance
    rays = uniform_rays(400_000, x=(-0.01, 0.01))
    full = 2.0 * np.radians(4.0) ** 2 * 9.9  # cm² sr Å, small angles
    expected = rays["p"].sum() / full
    value, counts = ps.brilliance(rays, ps.Window(0.005, 0.005, 1.0), np.array([0.1, 10.0]))
    assert value[0] == pytest.approx(expected, rel=0.03)
    assert counts[0] == pytest.approx(400_000 / 2 / 4, rel=0.03)


def test_the_transfer_needs_min_rays_in_both():
    value, reference = np.array([1.0, 1.0]), np.array([4.0, 4.0])
    rays, reference_rays = np.array([50.0, 5.0]), np.array([100.0, 100.0])
    transfer = ps.brilliance_transfer(value, rays, reference, reference_rays, min_rays=5)
    assert transfer.tolist() == [0.25, 0.25]
    transfer = ps.brilliance_transfer(value, rays, reference, reference_rays, min_rays=60)
    assert np.isnan(transfer).all()


def test_a_band_sums_the_bins_whose_centers_are_in_it():
    bins = np.array([0.0, 1.0, 2.0, 3.0])
    value, reference = np.array([1.0, 2.0, 4.0]), np.array([4.0, 4.0, 4.0])
    assert ps.band_transfer(value, reference, bins, (0.0, 2.0)) == pytest.approx(3.0 / 8.0)
    assert np.isnan(ps.band_transfer(value, np.zeros(3), bins, (0.0, 2.0)))
