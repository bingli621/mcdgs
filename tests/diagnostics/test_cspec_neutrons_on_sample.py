import numpy as np
import pytest

from mcdgs.diagnostics.cspec.neutrons_on_sample import CSPECToSample, subpulse_windows
from mcdgs.diagnostics.neutrons_on_sample import windows_table
from mcdgs.tof import build_tof


def test_the_instrument_ends_at_the_x_and_y_monitor_without_the_detector_panels():
    cspec = CSPECToSample()
    names = list(cspec.components)
    assert names[-2:] == ["sample_monitor", "sample_xy_monitor"]
    assert not any(name.startswith("detector_panel") for name in names)


def test_each_window_holds_the_elastic_arrival_and_its_energy_that_of_lambda():
    cspec = CSPECToSample()
    windows = subpulse_windows(cspec, build_tof(cspec, neutrons=100_000, seed=1).run())
    assert windows.sizes["subpulse"] >= 13
    at_sample = windows.coords["time_at_sample"].to(unit="ms").values
    assert np.all(np.diff(at_sample) > 0)
    opens = windows.coords["earliest_arrival"].to(unit="ms").values
    closes = windows.coords["latest_arrival"].to(unit="ms").values
    wavelength = windows.coords["incident_wavelength"].to(unit="Å").values
    elastic = at_sample + 1e3 * 3.5 / (3956.034 / wavelength)  # 3.5 m to the panels
    assert np.all((opens < elastic) & (elastic < closes))
    assert windows.coords["elastic_arrival"].to(unit="ms").values == pytest.approx(elastic)
    # the sub-pulse the choppers are phased to has lambda0
    assert np.min(np.abs(wavelength - cspec.params.lambda0)) < 0.01
    table = windows_table(windows).splitlines()
    assert table[0].startswith("| Sub-pulse | λi | Ei |")
    assert len(table) == 2 + windows.sizes["subpulse"]
