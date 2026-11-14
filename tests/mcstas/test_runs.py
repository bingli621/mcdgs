import json
from datetime import datetime, timedelta, timezone

import pytest

from mcdgs.mcstas.runs import Run, RunType

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
