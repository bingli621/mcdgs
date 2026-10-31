import re

import pytest

from mcdgs.mcstas.run_folders import (
    clean_runs,
    latest_run,
    list_runs,
    monitor_file,
    run_folder,
)


@pytest.mark.parametrize(
    "folder_name", ["", ".", "..", "../data", "a/b", "/tmp", ".hidden", "a b", 3]
)
def test_a_folder_name_is_one_safe_folder_name(tmp_path, folder_name):
    with pytest.raises(ValueError, match="folder_name"):
        run_folder(tmp_path / "out", folder_name, "run")
    assert not (tmp_path / "out").exists()


def test_runs_without_a_name_get_new_folders_by_time_and_kind(tmp_path):
    first = run_folder(tmp_path / "out", None, "scan_omega")
    assert (tmp_path / "out").is_dir() and not first.exists()  # mcrun creates the run folder
    first.mkdir()  # as the run fills it
    second = run_folder(tmp_path / "out", None, "scan_omega")
    assert re.fullmatch(r"\d{8}-\d{6}_scan_omega", first.name)
    assert second != first and second.name.startswith(first.name[:8])
    assert first.exists()  # a new run never cleans an earlier one


def test_a_named_run_replaces_its_folder(tmp_path):
    (tmp_path / "mine").mkdir()
    (tmp_path / "mine" / "mccode.h5").touch()
    assert run_folder(tmp_path, "mine", "run") == tmp_path / "mine"
    assert not (tmp_path / "mine").exists()


def test_latest_run_is_the_newest_finished_one_matching(tmp_path, finished_run):
    finished_run(tmp_path / "20261020-100000_scan_omega", 1)
    finished_run(tmp_path / "20261021-100000_scan_omega", 2)
    finished_run(tmp_path / "20261022-100000_run_to_sample", 3)
    (tmp_path / "20261023-100000_scan_omega").mkdir()  # unfinished
    assert list_runs(tmp_path) == [
        "20261022-100000_run_to_sample",
        "20261021-100000_scan_omega",
        "20261020-100000_scan_omega",
    ]
    assert latest_run(tmp_path) == tmp_path / "20261022-100000_run_to_sample"
    assert latest_run(tmp_path, "*_scan_omega") == tmp_path / "20261021-100000_scan_omega"
    with pytest.raises(FileNotFoundError):
        latest_run(tmp_path, "*_from_sample")


def test_no_output_folder_has_no_runs(tmp_path):
    assert list_runs(tmp_path / "missing") == []


def test_clean_keeps_the_newest_runs_the_beam_and_unfinished_folders(tmp_path, finished_run):
    beam = finished_run(tmp_path / "a_beam", 1, beam=True)
    runs = [finished_run(tmp_path / f"run{n}", 10 + n) for n in range(4)]
    unfinished = tmp_path / "unfinished"
    unfinished.mkdir()

    assert sorted(clean_runs(tmp_path, keep=2)) == runs[:2]
    assert sorted(tmp_path.iterdir()) == [beam, *runs[2:], unfinished]


def test_clean_never_removes_a_run_to_sample_run(tmp_path, finished_run):
    beams = [finished_run(tmp_path / f"beam{n}", n, beam=True) for n in range(3)]
    scans = [finished_run(tmp_path / f"scan{n}", 10 + n) for n in range(2)]
    assert clean_runs(tmp_path, keep=0) == scans[::-1]
    assert sorted(tmp_path.iterdir()) == beams


def test_a_run_from_the_sample_is_normalized_by_the_newest_beam_runs_monitor(
    tmp_path, finished_run
):
    finished_run(tmp_path / "1_run_to_sample", 1, beam=True)
    beam = finished_run(tmp_path / "2_run_to_sample", 2, beam=True)
    scan = finished_run(tmp_path / "3_scan_omega_from_sample", 3)
    assert monitor_file(scan) == beam / "mccode.h5"
    assert monitor_file(finished_run(tmp_path / "4_scan_omega", 4)) is None  # its own monitor


def test_a_run_from_the_sample_without_a_beam_run_raises(tmp_path, finished_run):
    with pytest.raises(FileNotFoundError, match="run_to_sample"):
        monitor_file(finished_run(tmp_path / "1_scan_omega_from_sample", 1))
