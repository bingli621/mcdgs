from pathlib import Path

import matplotlib.pyplot as plt

from mcdgs.diagnostics.cli import (
    markdown_table,
    output_folder,
    parse_out,
    remove_mcdisplay_folders,
    save_figure,
    save_text,
)


def test_out_is_read_from_the_command_line():
    assert parse_out("A module.\n\nMore.", ["--out", "docs/figures"]) == Path("docs/figures")


def test_without_out_there_is_none():
    assert parse_out(None, []) is None


def test_without_out_the_folder_is_in_the_results_of_the_workspace(tmp_path, monkeypatch):
    monkeypatch.setenv("MCDGS_WORKSPACE", str(tmp_path))
    folder = output_folder(None, "trex_choppers")
    assert folder == tmp_path / "results" / "trex_choppers"
    assert folder.is_dir()  # created with its parents


def test_out_is_created_with_its_parents(tmp_path, monkeypatch):
    monkeypatch.setenv("MCDGS_WORKSPACE", str(tmp_path / "workspace"))
    out = tmp_path / "a" / "b"
    assert output_folder(out, "trex_choppers") == out
    assert out.is_dir()
    assert not (tmp_path / "workspace").exists()  # the workspace is left alone


def test_the_mcdisplay_folders_of_a_build_folder_are_removed(tmp_path):
    left = [tmp_path / name for name in ("TRex_mcdisplay", "TRex_mcdisplay_0", "TRex_mcdisplay_12")]
    for folder in left:
        (folder / "mccode").mkdir(parents=True)
    kept = [tmp_path / "TRex.instr", tmp_path / "TRex_mcdisplay.txt", tmp_path / "TRex_db"]
    kept[0].touch()
    kept[1].touch()  # a file, not a folder of mcdisplay
    kept[2].mkdir()
    assert remove_mcdisplay_folders(tmp_path) == left
    assert sorted(tmp_path.iterdir()) == sorted(kept)
    assert remove_mcdisplay_folders(str(tmp_path)) == []


def test_a_saved_figure_is_closed_and_its_path_printed(tmp_path, capsys):
    fig, _ = plt.subplots()
    path = save_figure(fig, tmp_path / "figure.png", dpi=50)
    assert path.is_file()
    assert not plt.fignum_exists(fig.number)
    assert capsys.readouterr().out == f"Saved {path}\n"


def test_a_saved_text_has_its_path_printed(tmp_path, capsys):
    path = save_text("a table\n", tmp_path / "table.txt")
    assert path.read_text() == "a table\n"
    assert capsys.readouterr().out == f"Saved {path}\n"


def test_a_markdown_table_has_minus_signs_and_empty_cells():
    table = markdown_table(["Name", "Gap"], [["a", "-1.5 ms"], ["b", ""]])
    assert table.splitlines() == [
        "| Name | Gap |",
        "| ---- | --- |",
        "| a | \u22121.5 ms |",
        "| b | |",
    ]
