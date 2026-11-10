from pathlib import Path

import matplotlib.pyplot as plt

from mcdgs.diagnostics.cli import (
    markdown_table,
    output_folder,
    parse_out,
    save_figure,
    save_text,
)


def test_out_is_read_from_the_command_line():
    assert parse_out("A module.\n\nMore.", ["--out", "docs/figures"]) == Path("docs/figures")


def test_without_out_there_is_none():
    assert parse_out(None, []) is None


def test_the_default_folder_is_used_without_out(tmp_path):
    assert output_folder(None, tmp_path) == tmp_path


def test_out_is_created_with_its_parents(tmp_path):
    out = tmp_path / "a" / "b"
    assert output_folder(out, tmp_path) == out
    assert out.is_dir()


def test_the_default_folder_is_created_too(tmp_path):
    default = tmp_path / "mcstas_output" / "trex_choppers"
    assert output_folder(None, default) == default
    assert default.is_dir()


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
