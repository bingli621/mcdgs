from pathlib import Path

import numpy as np
import pytest
import scipp as sc

from mcdgs.mcstas import (
    numeric_parameters,
    read_events,
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


BINS = {"toa": sc.linspace("toa", 0.0, 0.015, 4, unit="s")}


def test_a_run_that_is_no_scan_has_no_scan_maps(tmp_path, finished_run):
    with pytest.raises(ValueError, match="not a scan"):
        read_scan_maps(finished_run(tmp_path / "run", 0), "sample_monitor", BINS)
