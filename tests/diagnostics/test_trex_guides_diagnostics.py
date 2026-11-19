import matplotlib

matplotlib.use("Agg")

from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np

from mcdgs.diagnostics.trex import guides_diagnostics
from mcdgs.instruments.trex.choppers import DISKS
from mcdgs.mcstas import McStasScriptDiagnostics


def fake_histogram(self, variable, bins, weighted=True):
    """Half the beam lost at each point, in every bin, with enough rays everywhere."""
    assert variable == "l"
    n = len(bins) - 1
    names = ["after source", "before bw1", "after guide_29", "after guide_50a", "before sample_pos"]
    return {name: np.full(n, (1e10 if weighted else 1e3) * 0.5**i) for i, name in enumerate(names)}


def test_main_saves_the_diagnostics_the_transmission_and_the_diagram(tmp_path, monkeypatch):
    runs = []

    def run(self, folder_name, **kwargs):
        runs.append((folder_name, kwargs, list(self._points)))
        return SimpleNamespace(folder=tmp_path / folder_name)  # as a record

    def run_analysis(self, folder_name, **kwargs):
        runs.append((folder_name, kwargs))
        return None, plt.figure()

    monkeypatch.setattr(McStasScriptDiagnostics, "run", run)
    monkeypatch.setattr(McStasScriptDiagnostics, "plot", lambda self: plt.figure())
    monkeypatch.setattr(McStasScriptDiagnostics, "run_analysis", run_analysis)
    monkeypatch.setattr(McStasScriptDiagnostics, "histogram", fake_histogram)
    paths = guides_diagnostics.main(out=tmp_path, ncount=1e3, analysis_ncount=1e4, seed=2)
    assert [p.name for p in paths] == [
        "trex_guides_diagnostics.png",
        "trex_guides_transmission.png",
        "trex_intensity.png",
    ]
    assert all(p.stat().st_size > 0 for p in paths)
    (name, kwargs, points), analysis = runs
    assert (name, kwargs) == ("trex_guides_diagnostics", {"ncount": 1e3, "seed": 2})
    assert [p["after"] or p["before"] for p in points] == [
        "source",
        "bw1",
        "guide_29",
        "guide_50a",
        "sample_pos",
    ]
    assert analysis == ("trex_guides_diagnostics_analysis", {"ncount": 1e4, "seed": 2})


def test_the_instrument_is_the_white_beam_with_every_chopper_parked():
    trex = guides_diagnostics.define_instrument()
    assert trex.params.lambda0 is None
    assert (trex.source.Lmin, trex.source.Lmax) == (0.1, 10.0)
    instr_param_dict = {p.name: p.value for p in trex.instr_params}
    for name, disk in DISKS.items():
        assert (instr_param_dict[f"nu_{name}"], instr_param_dict[f"delay_{name}"]) == (0, 0)
        assert instr_param_dict[f"park_{name}"] == disk.phase_from_mode("FF")


def test_the_transmission_has_no_value_where_too_few_rays_arrive():
    intensity = {"a": np.array([10.0, 10.0, 0.0]), "b": np.array([4.0, 1.0, 1.0])}
    rays = {"a": np.array([100.0, 100.0, 0.0]), "b": np.array([40.0, 5.0, 40.0])}
    result = guides_diagnostics.transmission(intensity, rays, "a", "b")
    assert result[0] == 0.4
    assert np.isnan(result[1:]).all()  # too few rays; nothing at the start


def test_the_transmission_figure_has_the_sections_and_the_whole_way():
    bins = guides_diagnostics.WAVELENGTH_BINS
    intensity = fake_histogram(None, "l", bins)
    rays = fake_histogram(None, "l", bins, weighted=False)
    fig = guides_diagnostics.plot_transmission(intensity, rays, bins, "title")
    sections, overall = fig.axes
    assert [t.get_text() for t in sections.get_legend().get_texts()] == [
        "source → bw1",
        "bw1 → end of the curve",
        "end of the curve → sample position",
    ]
    assert [t.get_text() for t in overall.get_legend().get_texts()] == [
        "source → sample position",
        "bw1 → sample position",
    ]
    assert overall.get_yscale() == "log" and fig.get_suptitle() == "title"
    plt.close(fig)
