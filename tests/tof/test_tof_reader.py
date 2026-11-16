import numpy as np
import pytest
import scipp as sc
import tof

from mcdgs.instruments.trex.model import TRex
from mcdgs.reduction.binning import monitor_bins
from mcdgs.tof import build_tof, read_tof_monitor


@pytest.fixture(scope="module")
def sample_monitor() -> tof.DetectorReading:
    return build_tof(TRex(), neutrons=20_000, seed=1).run()["sample_monitor"]


def _reached(monitor: tof.DetectorReading) -> sc.DataArray:
    neutrons = monitor.data["pulse", 0]
    return neutrons[~neutrons.masks["blocked_by_others"]]


def test_edges_in_another_unit_bin_as_in_the_unit_of_the_coordinate(sample_monitor):
    times_ms = np.linspace(0, 140, 71)  # tof keeps the time of arrival in µs
    counts = read_tof_monitor(
        sample_monitor, {"toa": sc.array(dims=["toa"], values=times_ms, unit="ms")}
    )
    expected, _ = np.histogram(
        _reached(sample_monitor).coords["toa"].to(unit="ms").values, times_ms
    )
    np.testing.assert_array_equal(counts, expected)


def test_only_the_neutrons_that_get_through_are_counted(sample_monitor):
    wavelengths = sc.linspace("wavelength", 0.0, 10.0, 51, unit="Å")
    counts = read_tof_monitor(sample_monitor, {"wavelength": wavelengths})
    blocked = sample_monitor.data["pulse", 0].masks["blocked_by_others"].values
    assert 0 < counts.sum() == np.count_nonzero(~blocked) < blocked.size


def test_there_is_one_axis_per_entry_of_bins_in_their_order(sample_monitor):
    bins = {
        "wavelength": sc.linspace("wavelength", 0.0, 10.0, 21, unit="Å"),
        "toa": sc.linspace("toa", 0.0, 0.14, 31, unit="s"),
    }
    assert read_tof_monitor(sample_monitor, bins).shape == (20, 30)


def test_the_sample_monitor_is_read_on_its_own_bins(sample_monitor):
    bins = monitor_bins(TRex().sample_monitor.bin_limits)
    shape = tuple(n_bins for _, _, n_bins in TRex().sample_monitor.bin_limits.values())
    assert read_tof_monitor(sample_monitor, bins).shape == shape  # time, then wavelength
