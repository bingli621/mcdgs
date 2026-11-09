from pathlib import Path

import h5py
import numpy as np
import pytest

from mcdgs.mcstas.nexus_io import merge_entries, numeric_parameters, read_events, read_pixels

SPLIT = Path(__file__).parents[1] / "data" / "split_run" / "split"  # from mcstas_to_sqe


def test_only_the_numeric_parameters_are_kept():
    assert numeric_parameters({"omega": "5", "file": "beam.mcpl", "x": None}) == {"omega": 5.0}


def test_single_run_files_become_one_entry_each_in_their_order(tmp_path):
    files = []
    for omega in (-5, 5):
        files.append(tmp_path / f"run_{omega}.h5")
        with h5py.File(files[-1], "w") as f:
            f.create_group("entry1").attrs["omega"] = omega

    with h5py.File(merge_entries(files, tmp_path / "scan.h5")) as scan:
        assert list(scan) == ["entry1", "entry2"]
        assert [scan[e].attrs["omega"] for e in scan] == [-5, 5]


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
