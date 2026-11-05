from pathlib import Path

from mcdgs.diagnostics.cli import output_folder, parse_out


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
