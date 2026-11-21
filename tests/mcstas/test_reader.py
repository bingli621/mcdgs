from pathlib import Path

import mcstasscript as ms
import numpy as np
import pytest
import scipp as sc
from mcstasscript.data.data import McStasData, McStasMetaData

from mcdgs.mcstas import (
    numeric_parameters,
    read_events,
    read_mcstas_map,
    read_pixels,
)

BINS = {
    "toa": sc.linspace("toa", 0.0, 15.0, 4, unit="ms"),
    "wavelength": sc.linspace("wavelength", 0.0, 10.0, 3, unit="Å"),
}


def _monitor(name: str, shape: tuple[int, int]) -> McStasData:
    """A McStas monitor as McStasScript has it: (wavelength, time)."""
    metadata = McStasMetaData()
    metadata.component_name = name
    data = McStasData(metadata)
    data.Intensity = np.arange(np.prod(shape), dtype=float).reshape(shape)
    data.Error = np.ones(shape)
    return data


@pytest.fixture
def loaded(tmp_path, finished_run, monkeypatch):
    """loaded(data) -> a single run whose folder ``ms.load_data`` reads as ``data``."""

    def make(data):
        run = finished_run(tmp_path / "run", 0)
        monkeypatch.setattr(
            ms, "load_data", lambda folder: data if folder == str(run.folder) else []
        )
        return run

    return make


def test_the_map_of_the_named_monitor_is_read_as_time_by_wavelength(loaded):
    data = [_monitor("other", (2, 3)), _monitor("sample_monitor", (2, 3))]
    data[0].Intensity += 100
    sample_map = read_mcstas_map("sample_monitor", loaded(data), BINS)
    assert sample_map.sizes == {"toa": 3, "wavelength": 2}  # time bins, wavelength bins
    np.testing.assert_array_equal(sample_map.values, [[0, 3], [1, 4], [2, 5]])


def test_the_map_is_in_n_per_s_with_its_errors_and_its_bin_edges(loaded):
    data = [_monitor("sample_monitor", (2, 3))]
    data[0].Error *= 2
    sample_map = read_mcstas_map("sample_monitor", loaded(data), BINS)
    assert sample_map.unit == "1/s"
    np.testing.assert_array_equal(sc.stddevs(sample_map).values, np.full((3, 2), 2.0))
    assert all(sc.identical(sample_map.coords[dim], BINS[dim]) for dim in BINS)
    assert sc.stddevs(sample_map.sum("toa")).values == pytest.approx([12**0.5] * 2)  # quadrature


def test_a_map_on_other_bins_is_refused(loaded):
    run = loaded([_monitor("sample_monitor", (3, 2))])
    with pytest.raises(ValueError, match=r"sample_monitor: McStas has a map of \(2, 3\) bins"):
        read_mcstas_map("sample_monitor", run, BINS)


def test_a_name_of_several_data_sets_is_refused(loaded):
    run = loaded([_monitor("sample_monitor", (2, 3)), _monitor("sample_monitor", (2, 3))])
    with pytest.raises(ValueError, match="sample_monitor: 2 data sets of that name"):
        read_mcstas_map("sample_monitor", run, BINS)


SPLIT = Path(__file__).parents[1] / "data" / "split_run" / "split"  # from mcstas_to_sqe


def test_only_the_numeric_parameters_are_kept():
    assert numeric_parameters({"omega": "5", "file": "beam.mcpl", "x": None}) == {"omega": 5.0}


def test_read_pixels_needs_the_beam_direction_for_a_file_without_a_source():
    with pytest.raises(ValueError, match="beam_direction"):
        read_pixels(SPLIT, "detector", "source", "sample_pos")


def test_read_pixels_puts_the_missing_source_upstream_along_the_beam_and_adds_no_coords():
    points = read_pixels(SPLIT, "detector", "source", "sample_pos", beam_direction=[0, 0, 2])
    assert [p["omega"] for p, _ in points] == [-5.0, 5.0]
    pixels = points[0][1]
    upstream = pixels.coords["sample_position"].value - pixels.coords["source_position"].value
    assert np.allclose(upstream, [0, 0, 1])  # 1 m, whatever the length of beam_direction
    assert "gamma" not in pixels.coords and "y" not in pixels.coords


def test_read_events_gives_each_neutron_its_time_in_ms_and_its_weight_in_n_per_s():
    binned = read_events(SPLIT, "detector", "sample_pos", "sample_pos")
    assert binned.dims == ("pixel_id",)
    assert binned.bins.unit == "1/s"
    assert binned.bins.coords["time_at_detector"].bins.unit == "ms"
    assert "t" not in binned.bins.coords
    # the same neutrons as read_pixels sums, of the first scan point
    (_, pixels), _ = read_pixels(
        SPLIT, "detector", "source", "sample_pos", beam_direction=[0, 0, 1]
    )
    assert np.allclose(binned.bins.sum().values, pixels.values)


def test_read_events_puts_the_missing_source_upstream_along_the_beam():
    with pytest.raises(ValueError, match="beam_direction"):
        read_events(SPLIT, "detector", "source", "sample_pos")
    binned = read_events(SPLIT, "detector", "source", "sample_pos", beam_direction=[0, 0, 3])
    upstream = binned.coords["sample_position"].value - binned.coords["source_position"].value
    assert np.allclose(upstream, [0, 0, 1])


def test_read_events_reads_each_scan_point_as_read_pixels_does():
    points = read_pixels(SPLIT, "detector", "source", "sample_pos", beam_direction=[0, 0, 1])
    for entry, (_, pixels) in enumerate(points, start=1):
        binned = read_events(
            SPLIT, "detector", "source", "sample_pos", beam_direction=[0, 0, 1], entry=entry
        )
        assert np.allclose(binned.bins.sum().values, pixels.values)
