import pytest

from mcdgs.reduction.binning import bin_edges


def test_the_edges_include_both_ends():
    edges = bin_edges("Qx", -1.0, 1.0, 0.1, unit="1/Å")
    assert edges.dims == ("Qx",)
    assert edges.unit == "1/Å"
    assert edges.sizes["Qx"] == 21
    assert edges.values[0] == pytest.approx(-1.0)
    assert edges.values[-1] == pytest.approx(1.0)


def test_without_a_step_one_bin_integrates():
    edges = bin_edges("l", -0.1, 0.1)
    assert edges.values.tolist() == pytest.approx([-0.1, 0.1])
    assert edges.unit == "dimensionless"


@pytest.mark.parametrize(
    "start, stop, step", [(1.0, 1.0, 0.1), (1.0, 0.0, 0.1), (0.0, 1.0, 0.0), (0.0, 1.0, 0.3)]
)
def test_a_bad_range_or_step_raises(start, stop, step):
    with pytest.raises(ValueError):
        bin_edges("h", start, stop, step)
