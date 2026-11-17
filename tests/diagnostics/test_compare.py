import matplotlib.pyplot as plt
import numpy as np
import pytest
import scipp as sc
import tof
from matplotlib.colors import LogNorm
from mcstasscript.data.data import McStasData, McStasMetaData

from mcdgs.diagnostics.compare import (
    intensity_label,
    mean_coordinate,
    plot_comparison,
    plot_map,
    plot_profile,
    read_comparison,
    read_scaled_comparisons,
)
from mcdgs.diagnostics.trex.source import TRexSource
from mcdgs.tof import build_tof, read_tof_map

BINS = {
    "toa": sc.linspace("toa", 0.0, 15.0, 31, unit="ms"),
    "wavelength": sc.linspace("wavelength", 0.0, 10.0, 21, unit="Å"),
}
NEUTRONS = 10_000


@pytest.fixture(scope="module")
def tof_result() -> tof.Result:
    return build_tof(TRexSource(), neutrons=NEUTRONS, seed=1).run()


def _mcstas_map(shape: tuple[int, int], name: str = "source_monitor") -> list[McStasData]:
    """A McStas monitor, source_monitor by default, as McStasScript has it: (wavelength, time)."""
    metadata = McStasMetaData()
    metadata.component_name = name
    data = McStasData(metadata)
    data.Intensity = np.arange(np.prod(shape), dtype=float).reshape(shape)
    data.Error = np.ones(shape)
    return [data]


def test_tof_is_scaled_to_the_total_of_mcstas(tof_result):
    comparison = read_comparison(
        "source_monitor", _mcstas_map((20, 30)), tof_result, BINS, 1e6, NEUTRONS
    )
    assert comparison.mcstas_intensity.shape == (30, 20)  # time, wavelength
    assert comparison.tof_intensity.sum() == pytest.approx(comparison.mcstas_intensity.sum())


def test_a_given_scale_is_kept(tof_result):
    comparison = read_comparison(
        "source_monitor", _mcstas_map((20, 30)), tof_result, BINS, 1e6, NEUTRONS, scale=2.0
    )
    counts = read_tof_map("source_monitor", tof_result, BINS).values
    np.testing.assert_allclose(comparison.tof_intensity, 2 * counts)
    np.testing.assert_allclose(comparison.tof_error, 2 * np.sqrt(counts))  # Poisson


def test_the_transmission_scales_tof_on_top_of_a_given_scale(tof_result):
    comparison = read_comparison(
        "source_monitor",
        _mcstas_map((20, 30)),
        tof_result,
        BINS,
        1e6,
        NEUTRONS,
        scale=2.0,
        transmission=0.25,
    )
    counts = read_tof_map("source_monitor", tof_result, BINS).values
    np.testing.assert_allclose(comparison.tof_intensity, 0.5 * counts)
    np.testing.assert_allclose(comparison.tof_error, 0.5 * np.sqrt(counts))
    assert (comparison.scale, comparison.transmission) == (2.0, 0.25)


def test_a_transmission_without_a_scale_is_refused(tof_result):
    with pytest.raises(ValueError, match="needs a given scale"):
        read_comparison(
            "source_monitor",
            _mcstas_map((20, 30)),
            tof_result,
            BINS,
            1e6,
            NEUTRONS,
            transmission=0.5,
        )


def test_a_map_on_other_bins_is_refused(tof_result):
    with pytest.raises(ValueError, match="source_monitor"):
        read_comparison("source_monitor", _mcstas_map((30, 20)), tof_result, BINS, 1e6, NEUTRONS)


def _two_monitors() -> list[McStasData]:
    """The same McStas map twice: as source_monitor and as a monitor further down the beam."""
    (first,) = _mcstas_map((20, 30))
    (second,) = _mcstas_map((20, 30), name="later_monitor")
    second.Intensity = first.Intensity / 4  # a quarter of the beam gets there
    return [first, second]


@pytest.fixture(scope="module")
def tof_two_monitors(tof_result) -> dict:
    """The tof reading of source_monitor, also as later_monitor."""
    return {
        "source_monitor": tof_result["source_monitor"],
        "later_monitor": tof_result["source_monitor"],
    }


def test_tof_is_scaled_at_the_first_monitor_and_the_transmission_tuned_at_the_second(
    tof_two_monitors,
):
    comparisons = read_scaled_comparisons(
        _two_monitors(),
        tof_two_monitors,
        {"source_monitor": BINS, "later_monitor": BINS},
        1e6,
        NEUTRONS,
    )
    first, later = comparisons["source_monitor"], comparisons["later_monitor"]
    assert list(comparisons) == ["source_monitor", "later_monitor"]
    assert first.transmission == 1.0
    assert later.scale == first.scale
    assert later.transmission == pytest.approx(0.25)
    assert later.tof_intensity.sum() == pytest.approx(later.mcstas_intensity.sum())


def test_a_given_transmission_is_kept_after_the_first_monitor(tof_two_monitors):
    comparisons = read_scaled_comparisons(
        _two_monitors(),
        tof_two_monitors,
        {"source_monitor": BINS, "later_monitor": BINS},
        1e6,
        NEUTRONS,
        transmission=0.5,
    )
    first, later = comparisons["source_monitor"], comparisons["later_monitor"]
    np.testing.assert_allclose(later.tof_intensity, 0.5 * first.tof_intensity)
    assert later.transmission == 0.5


def test_one_monitor_is_not_enough_to_scale(tof_result):
    with pytest.raises(ValueError, match="one after it"):
        read_scaled_comparisons(
            _mcstas_map((20, 30)), tof_result, {"source_monitor": BINS}, 1e6, NEUTRONS
        )


def test_the_figure_has_the_two_maps_the_spectrum_and_the_profile(tof_result):
    comparison = read_comparison(
        "source_monitor", _mcstas_map((20, 30)), tof_result, BINS, 1e6, NEUTRONS
    )
    fig = plot_comparison(comparison, "A title")
    assert fig.get_suptitle() == "A title"
    assert len(fig.axes) == 5  # four plots and the color bar
    plt.close(fig)


def test_the_maps_can_be_drawn_on_a_log_scale(tof_result):
    comparison = read_comparison(
        "source_monitor", _mcstas_map((20, 30)), tof_result, BINS, 1e6, NEUTRONS
    )
    fig = plot_comparison(comparison, "A title", log_decades=3, connect=True)
    mesh = fig.axes[0].collections[0]
    vmax = max(comparison.mcstas_intensity.max(), comparison.tof_intensity.max())
    assert isinstance(mesh.norm, LogNorm)
    assert mesh.norm.vmin == pytest.approx(vmax / 1e3)
    assert mesh.get_array().mask.any()  # the empty bins, in white
    assert mesh.cmap.name.startswith("viridis")
    plt.close(fig)


def test_the_spectrum_and_the_profile_can_be_on_a_log_axis(tof_result):
    comparison = read_comparison(
        "source_monitor", _mcstas_map((20, 30)), tof_result, BINS, 1e6, NEUTRONS
    )
    fig = plot_comparison(comparison, "A title", log_y=True)
    assert [ax.get_yscale() for ax in fig.axes[2:4]] == ["log", "log"]
    plt.close(fig)


def _map(values) -> sc.DataArray:
    """A map of 2 time bins by 3 wavelength bins, with errors of 1."""
    data = sc.array(
        dims=["toa", "wavelength"], values=np.asarray(values, float), variances=np.ones((2, 3))
    )
    edges = {
        "toa": sc.linspace("toa", 0.0, 0.002, 3, unit="s"),
        "wavelength": sc.linspace("wavelength", 1.0, 4.0, 4, unit="Å"),
    }
    return sc.DataArray(data, coords=edges)


def test_the_mean_of_a_coordinate_is_weighted_by_the_intensity():
    data = _map([[0, 0, 0], [1, 2, 1]])  # all in the second time bin, around 2.5 Å
    assert mean_coordinate(data, "toa").to(unit="ms").value == pytest.approx(1.5)
    assert mean_coordinate(data, "wavelength").value == pytest.approx(2.5)


def test_the_intensity_is_labeled_with_the_bin_widths_as_they_are_shown():
    data = _map(np.zeros((2, 3)))
    assert intensity_label(data) == "Intensity [n/s per 1 ms and 1 Å]"
    assert intensity_label(data.sum("toa")) == "Intensity [n/s per 1 Å]"


def test_a_map_is_drawn_with_its_first_dimension_along_x():
    fig, ax = plt.subplots()
    mesh = plot_map(ax, _map([[0, 1, 2], [3, 4, 5]]), vmax=10.0)
    assert (ax.get_xlabel(), ax.get_ylabel()) == ("Time [ms]", "Wavelength [Å]")
    assert mesh.get_array().shape == (3, 2)  # the rows along wavelength
    assert (mesh.norm.vmin, mesh.norm.vmax) == (0, 10.0)
    assert ax.get_xlim() == pytest.approx((0.0, 2.0))  # in ms
    plt.close(fig)


def test_a_profile_is_drawn_at_the_bin_centers_with_its_error_bars():
    fig, ax = plt.subplots()
    plot_profile(ax, _map([[0, 1, 2], [3, 4, 5]]).sum("toa"), "McStas")
    points, _, (bars,) = ax.containers[0]
    np.testing.assert_allclose(points.get_xdata(), [1.5, 2.5, 3.5])
    np.testing.assert_allclose(points.get_ydata(), [3, 5, 7])
    lower, upper = bars.get_segments()[0][:, 1]
    assert upper - lower == pytest.approx(2 * 2**0.5)  # two bins of error 1, in quadrature
    assert ax.get_xlabel() == "Wavelength [Å]"
    assert ax.get_ylabel() == "Intensity [n/s per 1 Å]"
    assert ax.get_legend_handles_labels()[1] == ["McStas"]
    plt.close(fig)
