import math

import numpy as np
import pytest
import scipp as sc

from mcdgs.diagnostics.cspec.elastic_resolution_incoherent_short import (
    FOCUS_MARGIN,
    define_instrument,
    detector_window,
)
from mcdgs.diagnostics.cspec.incident_energy import MC_MONITOR_TIME_RESOLUTION
from mcdgs.diagnostics.cspec.sample_runs import MONITOR, pixel_solid_angle
from mcdgs.instruments.cspec.detectors import (
    PANEL_ARM_ANGLES,
    PANEL_DISTANCE,
    PANEL_HEIGHT,
    PANEL_PIXELS,
    PANEL_WIDTH,
    PANEL_X_BINS,
    PANEL_Y_BINS,
    panel_angle,
)
from mcdgs.reduction.solid_angle import rectangle_solid_angle


def test_the_instrument_has_the_vanadium_the_panels_and_finer_bins_at_mc_monitor():
    cspec = define_instrument()
    names = list(cspec.components)
    assert names.index("sample_pos") < names.index("sample") < names.index("detector_panel_0")
    start, stop, n_bins = getattr(cspec, MONITOR).bin_limits["time"]
    assert (stop - start) / n_bins == pytest.approx(MC_MONITOR_TIME_RESOLUTION)
    assert type(cspec.source).__name__ == "Polygon_ESS_butterfly"


def test_the_window_covers_every_panel_with_its_margin():
    window = detector_window()
    center = math.degrees(math.atan2(window["target_x"], window["target_z"]))
    half_width = math.degrees(math.atan2(PANEL_WIDTH / 2, PANEL_DISTANCE))
    for i in range(len(PANEL_ARM_ANGLES)):
        offset = abs(panel_angle(i) - center) + half_width + FOCUS_MARGIN
        assert offset <= window["focus_aw"] / 2 + 1e-9
    elevation = math.degrees(math.atan2(PANEL_HEIGHT / 2, PANEL_DISTANCE))
    assert window["focus_ah"] / 2 == pytest.approx(elevation + FOCUS_MARGIN)


def test_the_pixels_of_a_panel_add_up_to_the_panel():
    index = 3
    phi = math.radians(panel_angle(index))
    center = PANEL_DISTANCE * np.array([math.sin(phi), 0.0, math.cos(phi)])
    across = np.array([math.cos(phi), 0.0, -math.sin(phi)])
    x = (np.arange(PANEL_X_BINS) + 0.5) / PANEL_X_BINS * PANEL_WIDTH - PANEL_WIDTH / 2
    y = (np.arange(PANEL_Y_BINS) + 0.5) / PANEL_Y_BINS * PANEL_HEIGHT - PANEL_HEIGHT / 2
    xx, yy = (a.ravel() for a in np.meshgrid(x, y))
    positions = center + xx[:, None] * across + yy[:, None] * np.array([0.0, 1.0, 0.0])
    sample_position = sc.vector([0.0, 0.0, 0.0], unit="m")
    binned = sc.DataArray(
        sc.zeros(dims=["pixel"], shape=[xx.size]),
        coords={
            "pixel_id": sc.array(dims=["pixel"], values=index * PANEL_PIXELS + np.arange(xx.size)),
            "position": sc.vectors(dims=["pixel"], values=positions, unit="m"),
            "sample_position": sample_position,
        },
    )
    pixels = pixel_solid_angle(binned)
    panel = rectangle_solid_angle(
        sc.vector(center, unit="m"),
        sample_position,
        sc.scalar(PANEL_WIDTH, unit="m"),
        sc.scalar(PANEL_HEIGHT, unit="m"),
        normal=sc.vector([math.sin(phi), 0.0, math.cos(phi)]),
    )
    assert pixels.sum().value == pytest.approx(panel.value, rel=1e-9)
