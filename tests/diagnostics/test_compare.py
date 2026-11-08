import matplotlib.pyplot as plt
import numpy as np
import pytest
import scipp as sc
import tof
from matplotlib.colors import LogNorm
from mcstasscript.data.data import McStasData, McStasMetaData

from mcdgs.diagnostics.compare import plot_comparison, read_comparison
from mcdgs.diagnostics.trex_01_source import TRexSource
from mcdgs.tof.builder import build_tof
from mcdgs.tof.reader import read_tof_monitor

BINS = {
    "toa": sc.linspace("toa", 0.0, 15.0, 31, unit="ms"),
    "wavelength": sc.linspace("wavelength", 0.0, 10.0, 21, unit="Å"),
}
NEUTRONS = 10_000


@pytest.fixture(scope="module")
def tof_result() -> tof.Result:
    return build_tof(TRexSource(), neutrons=NEUTRONS, seed=1).run()


def _mcstas_map(shape: tuple[int, int]) -> list[McStasData]:
    """A McStas monitor named source_monitor, as McStasScript has it: (wavelength, time)."""
    metadata = McStasMetaData()
    metadata.component_name = "source_monitor"
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
    counts = read_tof_monitor(tof_result["source_monitor"], BINS)
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
    counts = read_tof_monitor(tof_result["source_monitor"], BINS)
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
