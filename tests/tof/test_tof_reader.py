import numpy as np
import pytest
import scipp as sc
import tof

from mcdgs.instruments.trex.model import TRex
from mcdgs.reduction.binning import monitor_bins
from mcdgs.tof import build_tof, read_tof_map


@pytest.fixture(scope="module")
def tof_result() -> tof.Result:
    return build_tof(TRex(), neutrons=20_000, seed=1).run()


@pytest.fixture(scope="module")
def sample_monitor(tof_result) -> tof.DetectorReading:
    return tof_result["sample_monitor"]


def _reached(monitor: tof.DetectorReading) -> sc.DataArray:
    neutrons = monitor.data["pulse", 0]
    return neutrons[~neutrons.masks["blocked_by_others"]]


def test_edges_in_another_unit_bin_as_in_the_unit_of_the_coordinate(tof_result, sample_monitor):
    times_ms = np.linspace(0, 140, 71)  # tof keeps the time of arrival in µs
    edges = sc.array(dims=["toa"], values=times_ms, unit="ms")
    counts = read_tof_map("sample_monitor", tof_result, {"toa": edges})
    expected, _ = np.histogram(
        _reached(sample_monitor).coords["toa"].to(unit="ms").values, times_ms
    )
    np.testing.assert_array_equal(counts.values, expected)
    assert sc.identical(counts.coords["toa"], edges)  # in the unit given


def test_the_counts_have_their_poisson_variances(tof_result, sample_monitor):
    wavelengths = sc.linspace("wavelength", 0.0, 10.0, 51, unit="Å")
    counts = read_tof_map("sample_monitor", tof_result, {"wavelength": wavelengths})
    assert counts.unit == "counts"
    np.testing.assert_array_equal(counts.variances, counts.values)  # an error of √N


def test_only_the_neutrons_that_get_through_are_counted(tof_result, sample_monitor):
    wavelengths = sc.linspace("wavelength", 0.0, 10.0, 51, unit="Å")
    counts = read_tof_map("sample_monitor", tof_result, {"wavelength": wavelengths})
    blocked = sample_monitor.data["pulse", 0].masks["blocked_by_others"].values
    assert 0 < counts.values.sum() == np.count_nonzero(~blocked) < blocked.size


def test_there_is_one_axis_per_entry_of_bins_in_their_order(tof_result, sample_monitor):
    bins = {
        "wavelength": sc.linspace("wavelength", 0.0, 10.0, 21, unit="Å"),
        "toa": sc.linspace("toa", 0.0, 0.14, 31, unit="s"),
    }
    assert read_tof_map("sample_monitor", tof_result, bins).shape == (20, 30)


def test_the_sample_monitor_is_read_on_its_own_bins(tof_result, sample_monitor):
    bins = monitor_bins(TRex().sample_monitor.bin_limits)
    shape = tuple(n_bins for _, _, n_bins in TRex().sample_monitor.bin_limits.values())
    assert read_tof_map("sample_monitor", tof_result, bins).shape == shape  # time, then wavelength


def test_of_a_chopper_only_the_neutrons_it_lets_through_are_counted_by_default(tof_result):
    wavelengths = sc.linspace("wavelength", 0.0, 10.0, 51, unit="Å")
    chopper = tof_result["ps1"].data["pulse", 0]
    reached = ~chopper.masks["blocked_by_others"].values
    passed = reached & ~chopper.masks["blocked_by_me"].values
    transmitted = read_tof_map("ps1", tof_result, {"wavelength": wavelengths})
    arrived = read_tof_map("ps1", tof_result, {"wavelength": wavelengths}, transmitted=False)
    assert transmitted.values.sum() == np.count_nonzero(passed)
    assert arrived.values.sum() == np.count_nonzero(reached)
    assert 0 < np.count_nonzero(passed) < np.count_nonzero(reached)  # the chopper blocks some


def test_a_monitor_blocks_no_neutron(tof_result):
    wavelengths = sc.linspace("wavelength", 0.0, 10.0, 51, unit="Å")
    transmitted = read_tof_map("sample_monitor", tof_result, {"wavelength": wavelengths})
    arrived = read_tof_map(
        "sample_monitor", tof_result, {"wavelength": wavelengths}, transmitted=False
    )
    assert sc.identical(transmitted, arrived)
