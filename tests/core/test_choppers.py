from typing import Any

import pytest

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
