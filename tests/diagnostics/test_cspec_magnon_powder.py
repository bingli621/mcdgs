import numpy as np
import pytest

from mcdgs.diagnostics.cspec import magnon_powder as diagnostic
from mcdgs.diagnostics.cspec.elastic_resolution_incoherent_short import detector_window
from mcdgs.diagnostics.cspec.incident_energy import MC_MONITOR_TIME_RESOLUTION
from mcdgs.diagnostics.cspec.sample_runs import MONITOR
from mcdgs.mcstas import Workspace
from mcdgs.samples.inelastic.powder import magnon_band


def test_the_powder_reads_its_own_table_and_scatters_once_towards_the_panels(tmp_path):
    path = diagnostic.sqw_path(Workspace(tmp_path / "ws"))
    assert path.is_absolute() and path == (tmp_path / "ws" / "builds" / "cspec_magnon_powder.sqw")
    sample = diagnostic.sample(path)
    assert sample.Sqw_coh == str(path)
    assert (sample.name, sample.order, sample.p_interact) == ("sample", 1, diagnostic.P_INTERACT)
    assert sample.d_phi == pytest.approx(detector_window()["focus_ah"])
    assert sample.split == diagnostic.SAMPLE_SPLIT


def test_the_instrument_has_the_powder_the_panels_and_finer_bins_at_mc_monitor(tmp_path):
    cspec = diagnostic.CSPECMagnonPowder(sqw_file=tmp_path / "m.sqw")
    names = list(cspec.components)
    assert names.index(MONITOR) < names.index("sample_pos") < names.index("sample")
    assert names.index("sample") < names.index("detector_panel_0")
    assert cspec.sample.component_name == "Isotropic_Sqw"
    start, stop, n_bins = getattr(cspec, MONITOR).bin_limits["time"]
    assert (stop - start) / n_bins == pytest.approx(MC_MONITOR_TIME_RESOLUTION)


def test_the_magnon_runs_from_its_gap_to_its_top():
    low, high = magnon_band(
        np.array([0.0, np.sqrt(3) * np.pi / 4.0]),
        js=diagnostic.MAGNON_JS,
        gap=diagnostic.MAGNON_GAP,
    )
    assert low[0] == pytest.approx(1.0) and high[0] == pytest.approx(1.0)
    assert high[1] == pytest.approx(8.0, abs=0.05)  # the corner of the zone, |Q| = √3 π / a


def test_the_table_is_a_density_in_energy_on_both_sides_past_the_top(tmp_path):
    path = diagnostic.write_sqw(tmp_path / "m.sqw")
    lines = [line for line in path.read_text().splitlines() if not line.startswith("#")]
    q, energy = (np.array(lines[i].split(), dtype=float) for i in (0, 1))
    table = np.loadtxt(lines[2:])
    assert q[-1] == pytest.approx(diagnostic.SQW_Q_MAX)
    assert energy[1] - energy[0] == pytest.approx(diagnostic.SQW_E_STEP)
    assert energy[-1] > diagnostic.MAGNON_TOP and energy[0] == pytest.approx(-energy[-1])
    assert table.shape == (q.size, energy.size)
    loss = table[:, energy >= 0]  # one magnon on the loss side, and the same on the gain side
    assert loss.sum(axis=1) * diagnostic.SQW_E_STEP == pytest.approx(1, abs=2e-3)
