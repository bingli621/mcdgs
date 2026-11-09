import numpy as np
import pytest
import scipp as sc

from mcdgs.reduction.frame import beam_frame


def test_the_sample_is_the_origin_and_downstream_is_z():
    sample = sc.vector([1.0, 0.0, 10.0], unit="m")
    upstream = sc.vector([0.0, 0.0, 0.0], unit="m")
    beam = (sample - upstream) / sc.norm(sample - upstream)
    downstream = sample + sc.scalar(2.0, unit="m") * beam
    assert beam_frame(sample, sample, upstream).value == pytest.approx([0.0, 0.0, 0.0])
    assert beam_frame(downstream, sample, upstream).value == pytest.approx([0.0, 0.0, 2.0])


def test_up_stays_up_and_x_is_to_the_left_of_the_beam():
    sample = sc.vector([0.0, 0.0, 0.0], unit="m")
    upstream = sc.vector([0.0, 0.0, -5.0], unit="m")
    positions = sc.vectors(dims=["pixel"], values=[[0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], unit="m")
    local = beam_frame(positions, sample, upstream)
    assert np.allclose(local.values, [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0]])


def test_the_positions_keep_their_unit():
    local = beam_frame(
        sc.vector([0.0, 0.0, 300.0], unit="cm"),
        sc.vector([0.0, 0.0, 1.0], unit="m"),
        sc.vector([0.0, 0.0, 0.0], unit="m"),
    )
    assert local.unit == "cm"
    assert local.value == pytest.approx([0.0, 0.0, 200.0])
