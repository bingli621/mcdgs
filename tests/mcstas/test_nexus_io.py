from pathlib import Path

import h5py
import numpy as np
import pytest
import scipp as sc

from mcdgs.mcstas import (
    RunType,
    numeric_parameters,
    read_events,
    read_monitor_map,
    read_pixels,
    read_scan_maps,
)

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


BINS = {
    "toa": sc.linspace("toa", 0.0, 0.015, 4, unit="s"),
    "wavelength": sc.linspace("wavelength", 0.0, 10.0, 3, unit="Å"),
}


def _write_monitors(path, maps, name="sample_monitor"):
    """A file as McStas writes it, with one entry per (data, errors) map of the monitor `name`,
    after a monitor whose name ends the same, so that it is found by its whole name."""
    with h5py.File(path, "w") as file:
        for n, (data, errors) in enumerate(maps, start=1):
            components = file.create_group(f"entry{n}/instrument/components")
            other = components.create_group(f"0001_pre_{name}/output/pre_{name}_t_L")
            other["data"], other["errors"] = data + 100, errors
            time_map = components.create_group(f"0002_{name}/output/{name}_t_L")
            time_map["data"], time_map["errors"] = data, errors
    return path


def test_the_map_of_a_monitor_is_read_in_n_per_s_with_its_errors_on_its_bins(tmp_path):
    data = np.arange(6, dtype=float).reshape(3, 2)  # (time, wavelength), as McStas writes it
    file = _write_monitors(tmp_path / "mccode.h5", [(data, np.full((3, 2), 2.0))])
    sample_map = read_monitor_map(file, "sample_monitor", BINS)
    assert sample_map.sizes == {"toa": 3, "wavelength": 2} and sample_map.unit == "1/s"
    np.testing.assert_array_equal(sample_map.values, data)
    np.testing.assert_array_equal(sc.stddevs(sample_map).values, np.full((3, 2), 2.0))
    assert all(sc.identical(sample_map.coords[dim], BINS[dim]) for dim in BINS)


def test_a_map_on_other_bins_or_of_a_missing_monitor_is_refused(tmp_path):
    file = _write_monitors(tmp_path / "mccode.h5", [(np.ones((2, 3)), np.ones((2, 3)))])
    with pytest.raises(ValueError, match=r"sample_monitor: McStas has a map of \(2, 3\) bins"):
        read_monitor_map(file, "sample_monitor", BINS)
    with pytest.raises(ValueError, match="0 components 'beam_monitor'"):
        read_monitor_map(file, "beam_monitor", BINS)


def test_the_maps_of_a_scan_are_read_along_its_variable_in_its_order(tmp_path, finished_run):
    run = finished_run(tmp_path / "scan", 0, RunType.SCAN, variable="b_rot", values=list(range(11)))
    # 11 entries: by name, entry10 and entry11 would come before entry2
    maps = [(np.full((3, 2), float(n)), np.zeros((3, 2))) for n in range(11)]
    _write_monitors(run.data_file, maps)
    scan = read_scan_maps(run, "sample_monitor", BINS, unit="deg")
    assert scan.dims == ("b_rot", "toa", "wavelength")
    np.testing.assert_array_equal(scan.values[:, 0, 0], np.arange(11))
    assert sc.identical(scan.coords["b_rot"], sc.arange("b_rot", 11, unit="deg"))


def test_a_run_that_is_no_scan_has_no_scan_maps(tmp_path, finished_run):
    with pytest.raises(ValueError, match="not a scan"):
        read_scan_maps(finished_run(tmp_path / "run", 0), "sample_monitor", BINS)
