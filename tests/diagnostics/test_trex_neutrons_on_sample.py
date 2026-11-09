import numpy as np
import pytest

from mcdgs.diagnostics.trex_07_neutrons_on_sample import (
    TRexToSample,
    subpulse_windows,
    windows_table,
)
from mcdgs.tof.builder import build_tof


def test_each_window_holds_the_elastic_arrival_and_ends_before_the_next_opens():
    trex = TRexToSample()
    windows = subpulse_windows(trex, build_tof(trex, neutrons=100_000, seed=1).run())
    opens = windows.coords["earliest_arrival"].to(unit="ms").values
    closes = windows.coords["latest_arrival"].to(unit="ms").values
    at_sample = windows.coords["time_at_sample"].to(unit="ms").values
    speed = 3956.034 / windows.coords["incident_wavelength"].to(unit="Å").values  # m/s
    elastic = at_sample + 1e3 * 3.0 / speed  # 3 m to the banana detector
    assert np.all((opens < elastic) & (elastic < closes))
    assert windows.coords["elastic_arrival"].to(unit="ms").values == pytest.approx(elastic)
    assert np.all(closes[:-1] < opens[1:])
    table = windows_table(windows).splitlines()
    assert table[0].startswith("| Sub-pulse | λi | Ei |")
    assert len(table) == 2 + windows.sizes["subpulse"]
