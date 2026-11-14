from pathlib import Path

from mcdgs.mcstas.workspace import Workspace


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
