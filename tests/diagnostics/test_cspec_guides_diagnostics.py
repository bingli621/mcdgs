import matplotlib

matplotlib.use("Agg")

from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np

from mcdgs.core.components.source import L_MAX, L_MIN
from mcdgs.diagnostics.cspec import guides_diagnostics
from mcdgs.instruments.cspec.choppers import DISKS
from mcdgs.mcstas import McStasScriptDiagnostics

POINTS = [
    "after source",
    "before bw1",
    "after guide_06_13",
    "after guide_09_02",
    "before sample_pos",
]


def fake_histogram(self, variable, bins, weighted=True):
    """Half the beam lost at each point, in every bin, with enough rays everywhere."""
    assert variable == "l"
    n = len(bins) - 1
    return {name: np.full(n, (1e10 if weighted else 1e3) * 0.5**i) for i, name in enumerate(POINTS)}


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
        "cspec_guides_diagnostics.png",
        "cspec_guides_transmission.png",
        "cspec_intensity.png",
    ]
    assert all(p.stat().st_size > 0 for p in paths)
    (name, kwargs, points), analysis = runs
    assert (name, kwargs) == ("cspec_guides_diagnostics", {"ncount": 1e3, "seed": 2})
    assert [p["after"] or p["before"] for p in points] == [
        "source",
        "bw1",
        "guide_06_13",
        "guide_09_02",
        "sample_pos",
    ]
    assert analysis == ("cspec_guides_diagnostics_analysis", {"ncount": 1e4, "seed": 2})


def test_the_instrument_is_the_white_beam_with_every_chopper_parked():
    cspec = guides_diagnostics.define_instrument()
    assert cspec.params.lambda0 is None
    assert (cspec.source.Lmin, cspec.source.Lmax) == (L_MIN, L_MAX)
    instr_param_dict = {p.name: p.value for p in cspec.instr_params}
    for name, disk in DISKS.items():
        assert (instr_param_dict[f"nu_{name}"], instr_param_dict[f"delay_{name}"]) == (0, 0)
        assert instr_param_dict[f"park_{name}"] == disk.phase()


def test_the_points_are_the_ends_of_the_s_bender_and_of_the_straight_guide():
    components = list(guides_diagnostics.define_instrument().components)
    bender_end, straight_end = components.index("guide_06_13"), components.index("guide_09_02")
    assert components.index("bw1") < bender_end < components.index("bw3")
    assert components[straight_end + 1] == "sample_pos"
