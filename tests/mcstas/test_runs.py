import json
import re
from datetime import datetime, timedelta, timezone

import numpy as np
import pytest

from mcdgs.mcstas.runs import (
    Run,
    RunType,
    clean_runs,
    latest_run,
    list_runs,
    monitor_file,
    new_run_folder,
    recorded_instr_params,
)

TO_SAMPLE, FROM_SAMPLE = RunType.RUN_TO_SAMPLE, RunType.SCAN_FROM_SAMPLE

STARTED = datetime(2026, 11, 14, 10, 15, tzinfo=timezone(timedelta(hours=1)))


def _scan(folder):
    return Run(
        folder=folder,
        instrument="trex_la2ni7_scan",
        run_type=RunType.SCAN_FROM_SAMPLE,
        started=STARTED,
        instr_params={"lambda0": 2.5, "bender_in": 1, "beam": "beam.mcpl"},
        seed=1,
        variable="omega",
        values=[0.0, 1.0, 2.0],
        beam_run="trex_to_sample",
    )


def test_a_run_type_is_the_name_of_the_method_that_made_the_run():
    assert [run_type.value for run_type in RunType] == [
        "run",
        "scan",
        "run_to_sample",
        "run_from_sample",
        "scan_from_sample",
    ]


def test_a_saved_run_is_read_back_as_it_was(tmp_path):
    scan = _scan(tmp_path)
    assert scan.save() == tmp_path / "run.json"
    assert Run.load(tmp_path) == scan
    assert Run.load(str(tmp_path)) == scan


def test_a_run_without_the_optional_values_is_read_back_as_it_was(tmp_path):
    run = Run(tmp_path, "TRex", RunType.RUN, STARTED, {"omega": 5.0}, ncount=1e7)
    run.save()
    assert Run.load(tmp_path) == run


def test_the_folder_is_not_in_the_record_so_a_run_can_be_moved(tmp_path):
    old = tmp_path / "old"
    old.mkdir()
    _scan(old).save()
    assert "folder" not in json.loads((old / "run.json").read_text())
    new = old.rename(tmp_path / "new")
    moved = Run.load(new)
    assert moved.folder == new and moved.name == "new"
    assert moved.data_file == new / "mccode.h5"


def test_only_a_run_to_the_sample_has_a_beam_file(tmp_path):
    to_sample = Run(tmp_path, "TRex", RunType.RUN_TO_SAMPLE, STARTED, {}, ncount=1e8)
    assert to_sample.beam_file == tmp_path / "beam.mcpl.gz"
    assert _scan(tmp_path).beam_file is None


def test_a_folder_without_a_record_is_not_a_run(tmp_path):
    with pytest.raises(FileNotFoundError):
        Run.load(tmp_path)


@pytest.mark.parametrize(
    "folder_name", ["", ".", "..", "../data", "a/b", "/tmp", ".hidden", "a b", 3]
)
def test_a_folder_name_is_one_safe_folder_name(tmp_path, folder_name):
    with pytest.raises(ValueError, match="folder_name"):
        new_run_folder(tmp_path / "out", folder_name, "run")
    assert not (tmp_path / "out").exists()


def test_runs_without_a_name_get_new_folders_by_time_and_label(tmp_path):
    first = new_run_folder(tmp_path / "out", None, "scan_omega")
    assert (tmp_path / "out").is_dir() and not first.exists()  # mcrun creates the run folder
    first.mkdir()  # as the run fills it
    second = new_run_folder(tmp_path / "out", None, "scan_omega")
    assert re.fullmatch(r"\d{8}-\d{6}_scan_omega", first.name)
    assert second != first and second.name.startswith(first.name[:8])
    assert first.exists()  # a new run never cleans an earlier one


def test_a_named_run_replaces_its_folder(tmp_path, finished_run):
    finished_run(tmp_path / "mine", 1)
    assert new_run_folder(tmp_path, "mine", "run") == tmp_path / "mine"
    assert not (tmp_path / "mine").exists()


def test_the_finished_runs_are_listed_newest_first_by_their_start(tmp_path, finished_run):
    old = finished_run(tmp_path / "z_named_last", 1, RunType.SCAN)
    new = finished_run(tmp_path / "a_named_first", 2, RunType.SCAN)
    (tmp_path / "unfinished").mkdir()  # no record
    (tmp_path / "unfinished" / "mccode.h5").touch()
    assert list_runs(tmp_path) == [new, old]
    assert list_runs(str(tmp_path)) == [new, old]


def test_no_runs_folder_has_no_runs(tmp_path):
    assert list_runs(tmp_path / "missing") == []


def test_the_latest_run_is_the_newest_of_its_type_and_instrument(tmp_path, finished_run):
    finished_run(tmp_path / "a", 1, RunType.SCAN)
    scan = finished_run(tmp_path / "b", 2, RunType.SCAN)
    other = finished_run(tmp_path / "c", 3, RunType.SCAN, instrument="trex_extraction")
    beam = finished_run(tmp_path / "d", 4, TO_SAMPLE)
    assert latest_run(tmp_path) == beam
    assert latest_run(tmp_path, RunType.SCAN) == other
    assert latest_run(tmp_path, RunType.SCAN, "TRex") == scan
    assert latest_run(tmp_path, instrument="trex_extraction") == other
    with pytest.raises(FileNotFoundError, match="type scan_from_sample and instrument TRex"):
        latest_run(tmp_path, FROM_SAMPLE, "TRex")
    with pytest.raises(FileNotFoundError, match="no finished run in"):
        latest_run(tmp_path / "missing")


def test_a_run_from_the_sample_is_normalized_by_the_run_its_beam_is_from(tmp_path, finished_run):
    beam = finished_run(tmp_path / "beam", 1, TO_SAMPLE)
    finished_run(tmp_path / "newer_beam", 2, TO_SAMPLE)  # not the one it started from
    scan = finished_run(tmp_path / "scan", 3, FROM_SAMPLE, beam_run="beam")
    assert monitor_file(scan) == beam.folder / "mccode.h5"
    assert monitor_file(beam) == beam.folder / "mccode.h5"  # its own monitor


def test_a_run_from_the_sample_without_its_beam_run_raises(tmp_path, finished_run):
    with pytest.raises(FileNotFoundError, match="no record of the run to the sample"):
        monitor_file(finished_run(tmp_path / "unknown", 1, FROM_SAMPLE))
    with pytest.raises(FileNotFoundError):
        monitor_file(finished_run(tmp_path / "removed", 2, FROM_SAMPLE, beam_run="gone"))


def test_clean_keeps_the_newest_runs_the_beam_and_unfinished_folders(tmp_path, finished_run):
    beam = finished_run(tmp_path / "a_beam", 1, TO_SAMPLE)
    runs = [finished_run(tmp_path / f"run{n}", 10 + n).folder for n in range(4)]
    unfinished = tmp_path / "unfinished"
    unfinished.mkdir()
    (unfinished / "mccode.h5").touch()

    assert clean_runs(tmp_path, keep=2) == runs[1::-1]
    assert sorted(tmp_path.iterdir()) == [beam.folder, *runs[2:], unfinished]


def test_clean_without_a_type_never_removes_a_run_to_the_sample(tmp_path, finished_run):
    beams = [finished_run(tmp_path / f"beam{n}", n, TO_SAMPLE).folder for n in range(3)]
    scans = [finished_run(tmp_path / f"scan{n}", 10 + n, RunType.SCAN).folder for n in range(2)]
    assert clean_runs(tmp_path, keep=0) == scans[::-1]
    assert sorted(tmp_path.iterdir()) == beams


def test_clean_of_a_type_removes_only_the_older_runs_of_that_type(tmp_path, finished_run):
    beams = [finished_run(tmp_path / f"beam{n}", n, TO_SAMPLE).folder for n in range(3)]
    scans = [finished_run(tmp_path / f"scan{n}", 10 + n, RunType.SCAN).folder for n in range(3)]
    assert clean_runs(tmp_path, keep=1, run_type=TO_SAMPLE) == beams[1::-1]
    assert sorted(tmp_path.iterdir()) == [beams[2], *scans]
    assert clean_runs(tmp_path, keep=2, run_type=RunType.SCAN) == scans[:1]
    assert clean_runs(tmp_path, keep=5, run_type=RunType.RUN) == []
    assert sorted(tmp_path.iterdir()) == [beams[2], *scans[1:]]


def test_a_numpy_number_is_recorded_as_a_float():
    recorded = recorded_instr_params([], {"omega": np.float64(2.5), "n": np.int64(3)})
    assert recorded == {"omega": 2.5, "n": 3.0}
    assert all(type(value) is float for value in recorded.values())  # JSON takes them
