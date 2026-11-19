import json
from datetime import datetime, timedelta, timezone

import numpy as np
import pytest

from mcdgs.mcstas import Run, RunType, monitor_file, recorded_instr_params

TO_SAMPLE, FROM_SAMPLE = RunType.RUN_TO_SAMPLE, RunType.SCAN_FROM_SAMPLE

STARTED = datetime(2026, 11, 14, 10, 15, tzinfo=timezone(timedelta(hours=1)))


def _scan(folder):
    return Run(
        folder=folder,
        instrument="trex_la2ni7_scan",
        run_type=RunType.SCAN_FROM_SAMPLE,
        started=STARTED,
        instr_params={"lambda0": 2.5, "bender_in": 1, "beam": "beam.mcpl"},
        finished=STARTED + timedelta(minutes=90),
        seed=1,
        variable="omega",
        values=[0.0, 1.0, 2.0],
        beam_run="trex_to_sample",
    )


def test_a_run_type_is_the_name_of_the_method_that_made_the_run():
    assert [run_type.value for run_type in RunType] == [
        "run",
        "diagnose",
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


def test_the_duration_is_from_the_start_to_the_finish(tmp_path):
    assert _scan(tmp_path).duration == timedelta(minutes=90)
    assert Run(tmp_path, "TRex", RunType.RUN, STARTED, {}).duration is None  # not finished


def test_finish_saves_the_run_with_the_time_now_and_leaves_the_run_as_it_was(tmp_path):
    run = Run(tmp_path, "TRex", RunType.RUN, datetime.now().astimezone(), {"omega": 5.0})
    finished = run.finish()
    assert run.finished is None and not run.duration
    assert finished.finished.tzinfo is not None and finished.duration >= timedelta(0)
    assert Run.load(tmp_path) == finished


def test_a_record_from_before_the_finish_was_recorded_is_read_without_it(tmp_path):
    _scan(tmp_path).save()
    record = json.loads((tmp_path / "run.json").read_text())
    del record["finished"]
    (tmp_path / "run.json").write_text(json.dumps(record))
    old = Run.load(tmp_path)
    assert old.finished is None and old.duration is None
    assert old.started == STARTED


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


def test_a_numpy_number_is_recorded_as_a_float():
    recorded = recorded_instr_params([], {"omega": np.float64(2.5), "n": np.int64(3)})
    assert recorded == {"omega": 2.5, "n": 3.0}
    assert all(type(value) is float for value in recorded.values())  # JSON takes them
