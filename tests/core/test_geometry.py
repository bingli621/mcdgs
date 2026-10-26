import numpy as np
import pytest

from mcdgs.core.components import Arm
from mcdgs.core.geometry import place_absolute, rotation


def test_a_rotation_takes_the_lab_frame_to_the_rotated_one():
    assert rotation(0, 0, 0) == pytest.approx(np.eye(3))
    # rotated +90 deg about y: the frame's z is the lab's x, so lab x is z in the frame
    assert rotation(0, 90, 0) @ [1, 0, 0] == pytest.approx([0, 0, 1])


def test_place_absolute_follows_at_rotated_and_relative():
    components = {
        "a": Arm(at=[0, 0, 1]),
        "b": Arm(at=[0, 0, 0], rotated=[0, "omega", 0], relative="a"),
        "c": Arm(at=[0, 0, 2], relative="b"),  # 2 m along b's z
        "d": Arm(at=[0, 0, 1], rotated=[0, "-omega", 0], relative="c"),
    }
    placed = place_absolute(components, {"omega": 90})
    assert placed["c"][0] == pytest.approx([2, 0, 1])  # b's z is the lab's x
    assert placed["d"][0] == pytest.approx([3, 0, 1])
    assert placed["d"][1] == pytest.approx(np.eye(3))  # turned back
