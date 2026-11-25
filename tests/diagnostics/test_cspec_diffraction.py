import pytest

from mcdgs.diagnostics.cspec.diffraction import SAMPLE_SPLIT, define_instrument
from mcdgs.diagnostics.cspec.incident_energy import MC_MONITOR_TIME_RESOLUTION
from mcdgs.diagnostics.cspec.sample_runs import MONITOR


def test_the_instrument_has_the_crystal_the_panels_and_finer_bins_at_mc_monitor():
    cspec = define_instrument()
    names = list(cspec.components)
    assert names.index(MONITOR) < names.index("sample_pos") < names.index("sample")
    assert names.index("sample") < names.index("detector_panel_0")
    assert type(cspec.sample).__name__ == "Single_crystal"
    assert cspec.sample.split == SAMPLE_SPLIT
    start, stop, n_bins = getattr(cspec, MONITOR).bin_limits["time"]
    assert (stop - start) / n_bins == pytest.approx(MC_MONITOR_TIME_RESOLUTION)
    assert type(cspec.source).__name__ == "Polygon_ESS_butterfly"
