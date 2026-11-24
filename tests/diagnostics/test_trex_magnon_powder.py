import numpy as np
import pytest
import scipp as sc

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


def test_the_instrument_has_the_powder_on_the_sample_position(tmp_path):
    trex = diagnostic.TRexMagnonPowder(sqw_file=tmp_path / "m.sqw")
    assert trex.sample.component_name == "Isotropic_Sqw"
    assert trex.positions()["sample"] == trex.positions()["sample_pos"]
    assert diagnostic.define_instrument("RR").params.mode == "RR"


def test_the_maps_are_cropped_to_the_band_in_energy():
    edges = sc.linspace("energy_transfer", -40.0, 40.0, 801, unit="meV")
    q = sc.linspace("Q", 0.0, 5.0, 11, unit="1/Å")
    maps = sc.DataArray(
        sc.ones(dims=["subpulse", "energy_transfer", "Q"], shape=[2, 800, 10]),
        coords={"energy_transfer": edges, "Q": q},
    )
    cropped = diagnostic.crop(maps)
    energy = cropped.coords["energy_transfer"].values
    assert energy[0] == pytest.approx(diagnostic.ENERGY_RANGE[0])
    assert energy[-1] == pytest.approx(diagnostic.ENERGY_RANGE[1])
    assert cropped.sizes["Q"] == 10


def test_the_figure_has_a_panel_per_sub_pulse_with_the_band_over_it():
    n = 15
    edges = sc.linspace("energy_transfer", -10.0, 30.0, 401, unit="meV")
    q = sc.linspace("Q", 0.0, 6.0, 61, unit="1/Å")
    maps = sc.DataArray(
        sc.array(
            dims=["subpulse", "energy_transfer", "Q"],
            values=np.random.default_rng(1).random((n, 400, 60)),
        ),
        coords={
            "energy_transfer": edges,
            "Q": q,
            "incident_energy": sc.linspace("subpulse", 47.0, 11.0, n, unit="meV"),
        },
    )
    fig = diagnostic.plot_magnon_maps(maps, "test")
    panels = [ax for ax in fig.axes if ax.get_title().split(":")[0].isdigit()]
    assert len(panels) == n
    assert all(len(ax.lines) == 2 for ax in panels)  # the two edges of the band
