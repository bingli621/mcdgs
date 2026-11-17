import pytest

from mcdgs.diagnostics.drawing import component_range, remove_mcdisplay_folders, write_printout
from mcdgs.instruments.trex.model import TRex


def test_the_printout_lists_every_component(tmp_path):
    path = write_printout(TRex(), tmp_path / "trex_print.txt")
    assert path == tmp_path / "trex_print.txt"
    assert "sample_monitor" in path.read_text()


def test_the_components_between_two_leave_out_both():
    names = list(TRex().components)
    first, stop = component_range(TRex(), after="source", before="banana_detector")
    assert names[first - 1] == "source"
    assert names[stop] == "banana_detector"
    assert names[first:stop][0] == "source_monitor"
    assert names[first:stop][-1] == "beamstop"


def test_a_component_range_of_an_unknown_component_raises():
    with pytest.raises(ValueError, match="detector"):
        component_range(TRex(), after="source", before="detector")


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
