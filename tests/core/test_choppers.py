from typing import Any

import pytest

from mcdgs.core.components import Guide_gravity
from mcdgs.core.components.choppers import Disk


def disk(**changes: Any) -> Disk:
    values: dict[str, Any] = {
        "slits": ((0, 20), (125, 35)),
        "radius": 0.35,
        "beam_angle": 180,
        "zero_angle": 0,
        "slit_height": 0.095,
        "gap_arm": "arm",
        "gap_z": 100.0,
        "z": 0.05,
    }
    return Disk(**(values | changes))


def test_the_slit_edges_are_the_centers_plus_and_minus_half_the_widths():
    assert disk().slit_edges == (-10, 10, 107.5, 142.5)


def test_the_distance_is_the_arm_of_the_gap_plus_the_place_in_the_gap():
    assert disk().distance == pytest.approx(100.05)


def test_by_default_there_is_no_guide_in_front():
    assert disk().upstream_guide is None


def test_the_window_is_as_wide_as_the_exit_of_the_guide_in_front():
    guides = {"guide": Guide_gravity(w1=0.03, h1=0.1, w2=0.02, l=1)}
    assert disk(upstream_guide="guide").window_width(guides) == 0.02


def test_without_guides_the_window_is_unrestricted():
    assert disk(upstream_guide="guide").window_width(None) == 0


def test_without_a_guide_in_front_the_window_is_unrestricted():
    guides = {"guide": Guide_gravity(w1=0.03, h1=0.1, w2=0.02, l=1)}
    assert disk().window_width(guides) == 0
