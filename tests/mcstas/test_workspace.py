import re
from pathlib import Path

import pytest

from mcdgs.mcstas import RunType, Workspace

TO_SAMPLE, FROM_SAMPLE = RunType.RUN_TO_SAMPLE, RunType.SCAN_FROM_SAMPLE


def test_the_root_is_the_one_given(tmp_path, monkeypatch):
    monkeypatch.setenv("MCDGS_WORKSPACE", str(tmp_path / "variable"))
    assert Workspace.find(tmp_path / "given").root == tmp_path / "given"
    assert Workspace.find("~/given").root == Path.home() / "given"


def test_without_a_root_it_is_the_one_of_the_environment_variable(tmp_path, monkeypatch):
    monkeypatch.setenv("MCDGS_WORKSPACE", str(tmp_path / "variable"))
    assert Workspace.find().root == tmp_path / "variable"


def test_without_the_variable_the_root_is_in_the_current_folder(tmp_path, monkeypatch):
    monkeypatch.delenv("MCDGS_WORKSPACE", raising=False)
    monkeypatch.chdir(tmp_path)
    workspace = Workspace.find()
    assert workspace.root == tmp_path / "mcdgs_workspace"
    monkeypatch.chdir(tmp_path.parent)
    assert workspace.root == tmp_path / "mcdgs_workspace"  # fixed when it was found


def test_each_kind_of_file_has_its_folder_in_the_root(tmp_path):
    workspace = Workspace(tmp_path)
    assert workspace.builds_folder == tmp_path / "builds"
    assert workspace.runs_folder == tmp_path / "runs"
    assert workspace.results_folder == tmp_path / "results"
    assert workspace.build_folder("TRex") == tmp_path / "builds" / "TRex"
    assert workspace.run_folder("scan") == tmp_path / "runs" / "scan"
    assert workspace.result_folder("choppers") == tmp_path / "results" / "choppers"
    assert not any(tmp_path.iterdir())  # only paths: nothing is created


def test_a_root_given_as_a_string_is_a_path(tmp_path):
    assert Workspace(str(tmp_path)).root == tmp_path  # pyright: ignore[reportArgumentType]


def test_cleaning_the_builds_removes_every_build_folder_and_nothing_else(tmp_path):
    workspace = Workspace(tmp_path)
    assert workspace.clean_builds() == []  # no builds folder yet
    builds = [workspace.build_folder(name) for name in ("TRex", "trex_extraction")]
    for build in builds:
        (build / "TRex_mcdisplay_0").mkdir(parents=True)
    run = workspace.run_folder("scan")
    run.mkdir(parents=True)
    assert workspace.clean_builds() == builds
    assert workspace.builds_folder.is_dir() and not any(workspace.builds_folder.iterdir())
    assert run.is_dir()


@pytest.mark.parametrize(
    "folder_name", ["", ".", "..", "../data", "a/b", "/tmp", ".hidden", "a b", 3]
)
def test_a_folder_name_is_one_safe_folder_name(tmp_path, folder_name):
    with pytest.raises(ValueError, match="folder_name"):
        Workspace(tmp_path).new_run_folder(folder_name, "run")
    assert not any(tmp_path.iterdir())


def test_runs_without_a_name_get_new_folders_by_time_and_label(tmp_path):
    workspace = Workspace(tmp_path)
    first = workspace.new_run_folder(None, "scan_omega")
    assert first.parent == workspace.runs_folder
    assert workspace.runs_folder.is_dir() and not first.exists()  # mcrun creates the run folder
    first.mkdir()  # as the run fills it
    second = workspace.new_run_folder(None, "scan_omega")
    assert re.fullmatch(r"\d{8}-\d{6}_scan_omega", first.name)
    assert second != first and second.name.startswith(first.name[:8])
    assert first.exists()  # a new run never cleans an earlier one


def test_a_named_run_replaces_its_folder(tmp_path, finished_run):
    workspace = Workspace(tmp_path)
    finished_run(workspace.run_folder("mine"), 1)
    assert workspace.new_run_folder("mine", "run") == workspace.run_folder("mine")
    assert not workspace.run_folder("mine").exists()


def test_the_finished_runs_are_listed_newest_first_by_their_start(tmp_path, finished_run):
    workspace = Workspace(tmp_path)
    old = finished_run(workspace.run_folder("z_named_last"), 1, RunType.SCAN)
    new = finished_run(workspace.run_folder("a_named_first"), 2, RunType.SCAN)
    unfinished = workspace.run_folder("unfinished")  # no record
    unfinished.mkdir()
    (unfinished / "mccode.h5").touch()
    assert workspace.list_runs() == [new, old]


def test_no_runs_folder_has_no_runs(tmp_path):
    assert Workspace(tmp_path / "missing").list_runs() == []


def test_the_latest_run_is_the_newest_of_its_type_and_instrument(tmp_path, finished_run):
    workspace = Workspace(tmp_path)
    finished_run(workspace.run_folder("a"), 1, RunType.SCAN)
    scan = finished_run(workspace.run_folder("b"), 2, RunType.SCAN)
    other = finished_run(workspace.run_folder("c"), 3, RunType.SCAN, instrument="trex_extraction")
    beam = finished_run(workspace.run_folder("d"), 4, TO_SAMPLE)
    assert workspace.latest_run() == beam
    assert workspace.latest_run(RunType.SCAN) == other
    assert workspace.latest_run(RunType.SCAN, "TRex") == scan
    assert workspace.latest_run(instrument="trex_extraction") == other
    with pytest.raises(FileNotFoundError, match="type scan_from_sample and instrument TRex"):
        workspace.latest_run(FROM_SAMPLE, "TRex")
    with pytest.raises(FileNotFoundError, match="no finished run in"):
        Workspace(tmp_path / "missing").latest_run()


def test_clean_keeps_the_newest_runs_the_beam_and_unfinished_folders(tmp_path, finished_run):
    workspace = Workspace(tmp_path)
    beam = finished_run(workspace.run_folder("a_beam"), 1, TO_SAMPLE)
    runs = [finished_run(workspace.run_folder(f"run{n}"), 10 + n).folder for n in range(4)]
    unfinished = workspace.run_folder("unfinished")
    unfinished.mkdir()
    (unfinished / "mccode.h5").touch()

    assert workspace.clean_runs(keep=2) == runs[1::-1]
    assert sorted(workspace.runs_folder.iterdir()) == [beam.folder, *runs[2:], unfinished]


def test_clean_without_a_type_never_removes_a_run_to_the_sample(tmp_path, finished_run):
    workspace = Workspace(tmp_path)
    beams = [finished_run(workspace.run_folder(f"beam{n}"), n, TO_SAMPLE).folder for n in range(3)]
    scans = [
        finished_run(workspace.run_folder(f"scan{n}"), 10 + n, RunType.SCAN).folder
        for n in range(2)
    ]
    assert workspace.clean_runs(keep=0) == scans[::-1]
    assert sorted(workspace.runs_folder.iterdir()) == beams


def test_clean_of_a_type_removes_only_the_older_runs_of_that_type(tmp_path, finished_run):
    workspace = Workspace(tmp_path)
    beams = [finished_run(workspace.run_folder(f"beam{n}"), n, TO_SAMPLE).folder for n in range(3)]
    scans = [
        finished_run(workspace.run_folder(f"scan{n}"), 10 + n, RunType.SCAN).folder
        for n in range(3)
    ]
    assert workspace.clean_runs(keep=1, run_type=TO_SAMPLE) == beams[1::-1]
    assert sorted(workspace.runs_folder.iterdir()) == [beams[2], *scans]
    assert workspace.clean_runs(keep=2, run_type=RunType.SCAN) == scans[:1]
    assert workspace.clean_runs(keep=5, run_type=RunType.RUN) == []
    assert sorted(workspace.runs_folder.iterdir()) == [beams[2], *scans[1:]]
