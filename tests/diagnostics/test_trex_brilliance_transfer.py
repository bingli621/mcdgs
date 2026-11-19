import matplotlib

matplotlib.use("Agg")

from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

from mcdgs.diagnostics.trex import brilliance_transfer as bt
from mcdgs.instruments.trex.sources import source
from mcdgs.mcstas import McStasScriptDiagnostics


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
    window = bt.Window(half_width=0.01, half_height=0.005, half_divergence=0.5, x_center=0.05)
    rays = {
        "x": np.array([0.05, 0.065, 0.05, 0.05]),
        "y": np.zeros(4),
        "dx": np.array([0.0, 0.0, 0.6, 0.0]),
        "dy": np.array([0.0, 0.0, 0.0, -0.6]),
    }
    assert window.contains(rays).tolist() == [True, False, False, False]
    assert window.area == pytest.approx(2.0)  # cm²: 2 × 1 cm


def test_the_brilliance_of_an_even_beam_is_its_intensity_over_its_phase_space():
    # 2 cm × 1 cm, ±2° both ways, 9.9 Å: every window inside it has the same brilliance
    rays = uniform_rays(400_000, x=(-0.01, 0.01))
    full = 2.0 * np.radians(4.0) ** 2 * 9.9  # cm² sr Å, small angles
    expected = rays["p"].sum() / full
    value, counts = bt.brilliance(rays, bt.Window(0.005, 0.005, 1.0), np.array([0.1, 10.0]))
    assert value[0] == pytest.approx(expected, rel=0.03)
    assert counts[0] == pytest.approx(400_000 / 2 / 4, rel=0.03)


def test_each_face_is_read_over_its_span_on_the_plane():
    rays = uniform_rays(200_000)
    rays["p"] = np.where(rays["x"] < 0, 2.0, 1.0)  # the thermal side twice as bright
    faces = bt.moderator_brilliance(rays, np.array([0.1, 10.0]))
    assert list(faces) == ["thermal", "cold"]
    thermal, cold = faces["thermal"][0][0], faces["cold"][0][0]
    assert thermal / cold == pytest.approx(2.0, rel=0.05)


def test_the_reference_is_the_brighter_face_in_each_bin():
    faces = {
        "thermal": (np.array([3.0, 1.0]), np.array([30.0, 10.0])),
        "cold": (np.array([1.0, 2.0]), np.array([20.0, 40.0])),
    }
    value, counts = bt.reference_brilliance(faces)
    assert value.tolist() == [3.0, 2.0]
    assert counts.tolist() == [30.0, 40.0]


def test_a_band_sums_the_bins_whose_centers_are_in_it():
    bins = np.array([0.0, 1.0, 2.0, 3.0])
    value, reference = np.array([1.0, 2.0, 4.0]), np.array([4.0, 4.0, 4.0])
    assert bt.band_transfer(value, reference, bins, (0.0, 2.0)) == pytest.approx(3.0 / 8.0)
    assert np.isnan(bt.band_transfer(value, np.zeros(3), bins, (0.0, 2.0)))


def test_the_moderator_view_is_the_source_of_trex_aimed_at_its_plane():
    view = bt.ModeratorView()
    assert list(view.components) == ["origin", "source", "plane"]
    assert view.source == replace(
        source(),
        dist=bt.REFERENCE_DISTANCE,
        focus_xw=bt.REFERENCE_WIDTH,
        focus_yh=bt.REFERENCE_HEIGHT,
    )
    assert view.plane.at == [0, 0, bt.REFERENCE_DISTANCE]
    assert view.plane.relative == "source"
    # both faces, seen through the reference divergence, lie on the rectangle the source aims at
    reach = bt.REFERENCE_DISTANCE * np.tan(np.radians(bt.REFERENCE_HALF_DIVERGENCE))
    assert -bt.REFERENCE_WIDTH / 2 <= bt.THERMAL_FACE[0] < bt.THERMAL_FACE[1] < -reach
    assert reach < bt.COLD_FACE[0] < bt.COLD_FACE[1] <= bt.REFERENCE_WIDTH / 2


def test_main_runs_the_moderator_then_trex_and_saves_the_figure(tmp_path, monkeypatch):
    runs = []

    def run(self, folder_name, **kwargs):
        runs.append((self.simulation.instrument.name, folder_name, kwargs, list(self._points)))
        return SimpleNamespace(folder=tmp_path / folder_name)

    def rays(self, *variables):
        assert variables == bt.RAY_VARIABLES
        labels = [f"before {point['before']}" for point in self._points]
        return {label: uniform_rays(50_000, seed=i) for i, label in enumerate(labels)}

    monkeypatch.setattr(McStasScriptDiagnostics, "run", run)
    monkeypatch.setattr(McStasScriptDiagnostics, "rays", rays)
    path = bt.main(out=tmp_path, ncount=1e3, reference_ncount=1e2, seed=2)
    assert path == tmp_path / "trex_brilliance_transfer.png"
    assert path.stat().st_size > 0
    (view, view_folder, view_kwargs, view_points), (trex, folder, kwargs, points) = runs
    assert (view, view_folder, view_kwargs) == (
        "trex_moderator_view",
        "trex_moderator_brilliance",
        {"ncount": 1e2, "seed": 2},
    )
    assert [p["before"] for p in view_points] == ["plane"]
    assert (trex, folder, kwargs) == (
        "TRex",
        "trex_brilliance_transfer",
        {"ncount": 1e3, "seed": 2},
    )
    assert [p["before"] for p in points] == ["beam_monitor", "sample_pos"]


def test_the_instrument_is_that_of_guides_diagnostics():
    from mcdgs.diagnostics.trex import guides_diagnostics

    assert bt.define_instrument() == guides_diagnostics.define_instrument()
