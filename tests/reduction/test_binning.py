import pytest
import scipp as sc

from mcdgs.reduction.binning import bin_edges, monitor_bins


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


def test_the_bins_of_a_monitor_are_those_of_its_limits_in_their_order():
    bins = monitor_bins({"lambda": (0.0, 10.0, 20), "time": (0.0, 0.015, 30)})
    assert list(bins) == ["wavelength", "toa"]
    assert sc.identical(bins["toa"], sc.linspace("toa", 0.0, 0.015, 31, unit="s"))
    assert sc.identical(bins["wavelength"], sc.linspace("wavelength", 0.0, 10.0, 21, unit="Å"))


def test_an_axis_that_is_no_coordinate_of_the_neutrons_is_left_out():
    bins = monitor_bins({"x": (-0.05, 0.05, 10), "lambda": (1.0, 5.0, 4)})
    assert list(bins) == ["wavelength"]
    assert monitor_bins({}) == {}
