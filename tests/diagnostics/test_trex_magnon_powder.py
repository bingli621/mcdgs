import pytest

from mcdgs.diagnostics.trex import magnon_powder as diagnostic
from mcdgs.diagnostics.trex.solid_angles import detector_window
from mcdgs.mcstas import Workspace


def test_the_powder_reads_the_table_of_the_workspace_and_scatters_once_towards_the_detector(
    tmp_path,
):
    path = diagnostic.sqw_path(Workspace(tmp_path / "ws"))
    assert path.is_absolute() and path == (tmp_path / "ws" / "builds" / "magnon_powder.sqw")
    sample = diagnostic.sample(path)
    assert sample.Sqw_coh == str(path)
    assert (sample.name, sample.order, sample.p_interact) == ("sample", 1, diagnostic.P_INTERACT)
    assert sample.d_phi == pytest.approx(detector_window()["focus_ah"])
    assert sample.split == diagnostic.SAMPLE_SPLIT


def test_the_instrument_has_the_powder_on_the_sample_position(tmp_path):
    trex = diagnostic.TRexMagnonPowder(sqw_file=tmp_path / "m.sqw")
    assert trex.sample.component_name == "Isotropic_Sqw"
    assert trex.positions()["sample"] == trex.positions()["sample_pos"]
    assert diagnostic.define_instrument("RR").params.mode == "RR"
