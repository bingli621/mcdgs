import pytest

from mcdgs.diagnostics.cspec import diffraction
from mcdgs.diagnostics.cspec.elastic_resolution_incoherent_short import detector_window
from mcdgs.diagnostics.cspec.incident_energy import MC_MONITOR_TIME_RESOLUTION
from mcdgs.diagnostics.cspec.powder_diffraction import SAMPLE_SPLIT, define_instrument
from mcdgs.diagnostics.cspec.sample_runs import MONITOR


def test_the_instrument_has_the_powder_the_panels_and_finer_bins_at_mc_monitor():
    cspec = define_instrument()
    names = list(cspec.components)
    assert names.index(MONITOR) < names.index("sample_pos") < names.index("sample")
    assert names.index("sample") < names.index("detector_panel_0")
    assert type(cspec.sample).__name__ == "PowderN"
    assert cspec.sample.reflections == "Al2O3_sapphire.laz"
    assert cspec.sample.split == SAMPLE_SPLIT
    start, stop, n_bins = getattr(cspec, MONITOR).bin_limits["time"]
    assert (stop - start) / n_bins == pytest.approx(MC_MONITOR_TIME_RESOLUTION)
    assert type(cspec.source).__name__ == "Polygon_ESS_butterfly"


def test_the_powder_spreads_its_cones_over_the_height_of_the_panels():
    assert define_instrument().sample.d_phi == pytest.approx(detector_window()["focus_ah"])


def test_the_beam_of_the_crystal_is_that_of_the_powder_up_to_the_sample():
    powder, crystal = define_instrument(), diffraction.define_instrument()
    names = list(powder.components)
    before = names[: names.index("sample_pos") + 1]
    assert before == list(crystal.components)[: len(before)]
    assert [powder.components[n] for n in before] == [crystal.components[n] for n in before]
    assert powder.params == crystal.params
