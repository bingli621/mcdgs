import matplotlib

matplotlib.use("Agg")

import pytest
from matplotlib.figure import Figure

from mcdgs.instruments.trex.model import TRex
from mcdgs.tof import build_tof, time_distance_diagram
from mcdgs.tof.time_distance import ZOOM_MARGIN, wavelength_range, zoom_between


def test_the_title_is_the_running_conditions():
    trex = TRex()
    trex.update_params(mode="RR", rrm=12)
    figure = time_distance_diagram(trex, 10_000, seed=1)
    assert isinstance(figure, Figure)
    assert figure.axes[0].get_title() == trex.running_conditions()


def test_draws_the_pulses_given():
    figure = time_distance_diagram(TRex(), 10_000, seed=1, n_pulses=2, blocked_rays=0)
    labels = [text.get_text() for text in figure.axes[0].texts]
    assert labels.count("Pulse") == 2


def test_between_zooms_in_on_the_flight_from_one_component_to_the_next():
    trex = TRex()
    figure = time_distance_diagram(
        trex, 100_000, seed=1, blocked_rays=0, between=("pre_sample_monitor", "detector")
    )
    model = build_tof(trex, neutrons=1000, seed=1)
    near = model.components["pre_sample_monitor"].distance.value
    far = model.components["detector"].distance.value
    low, high = figure.axes[0].get_ylim()
    assert low < near < low + 0.1 * (far - near) and high - 0.1 * (far - near) < far < high
    assert all(text.get_clip_on() for text in figure.axes[0].texts)


def test_zoom_between_spans_the_neutrons_that_reach_the_second():
    result = build_tof(TRex(), neutrons=100_000, seed=1).run()
    limits = zoom_between(result, "pre_sample_monitor", "detector")
    reached = result["detector"].data
    through = ~reached.masks["blocked_by_others"].values
    start = result["pre_sample_monitor"].data.coords["toa"].values[through].min()
    end = reached.coords["toa"].values[through].max()
    low, high = limits["xlim"]
    assert low < start and end < high
    assert high - low == pytest.approx((end - start) * (1 + 2 * ZOOM_MARGIN))


@pytest.mark.parametrize(
    ("first", "second", "match"),
    [
        ("nowhere", "detector", "not a component"),
        ("detector", "pre_sample_monitor", "not nearer the source"),
    ],
)
def test_zoom_between_raises_for_names_it_cannot_zoom_on(first, second, match):
    result = build_tof(TRex(), neutrons=1000, seed=1).run()
    with pytest.raises(ValueError, match=match):
        zoom_between(result, first, second)


def test_color_by_spans_the_color_scale_over_the_wavelengths_at_a_component():
    trex = TRex()
    figure = time_distance_diagram(
        trex, 100_000, seed=1, blocked_rays=0, color_by="pre_sample_monitor"
    )
    result = build_tof(trex, neutrons=100_000, seed=1).run()
    low, high = wavelength_range(result, "pre_sample_monitor", "detector")
    assert figure.axes[1].get_ylim() == pytest.approx((low, high))  # the colorbar
    raised = time_distance_diagram(
        trex, 100_000, seed=1, blocked_rays=0, color_by="pre_sample_monitor", color_headroom=1.0
    )
    assert raised.axes[1].get_ylim() == pytest.approx((low, high + 1.0))


def test_wavelength_range_raises_for_an_unknown_component():
    result = build_tof(TRex(), neutrons=1000, seed=1).run()
    with pytest.raises(ValueError, match="not a component"):
        wavelength_range(result, "nowhere", "detector")


def test_a_title_goes_above_the_running_conditions():
    trex = TRex()
    figure = time_distance_diagram(trex, 10_000, seed=1, title="Zoomed")
    assert figure.axes[0].get_title() == f"Zoomed\n{trex.running_conditions()}"


def test_a_percentile_ends_the_time_axis_before_the_slowest_neutrons():
    result = build_tof(TRex(), neutrons=100_000, seed=1).run()
    full = zoom_between(result, "pre_sample_monitor", "detector")
    shorter = zoom_between(result, "pre_sample_monitor", "detector", percentile=50)
    assert shorter["xlim"][1] < full["xlim"][1]
    assert shorter["ylim"] == full["ylim"]
