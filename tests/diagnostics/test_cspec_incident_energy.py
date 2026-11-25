import numpy as np
import pytest
import scipp as sc

from mcdgs.diagnostics import compare
from mcdgs.diagnostics.cspec.incident_energy import (
    MC_MONITOR_TIME_RESOLUTION,
    MONITOR,
    define_instrument,
    subpulse_flights,
)
from mcdgs.reduction.binning import monitor_bins
from mcdgs.tof import build_tof, read_tof_maps


def test_the_instrument_ends_at_mc_monitor_with_finer_time_bins_and_the_masked_source():
    cspec = define_instrument()
    assert list(cspec.components)[-1] == MONITOR
    start, stop, n_bins = cspec.mc_monitor.bin_limits["time"]
    assert (stop - start) / n_bins == pytest.approx(MC_MONITOR_TIME_RESOLUTION)
    assert type(cspec.source).__name__ == "Polygon_ESS_butterfly"


def test_the_flight_from_the_source_gives_the_energy_of_the_mean_wavelength():
    cspec = define_instrument()
    bins = {MONITOR: monitor_bins(cspec.mc_monitor.bin_limits)}
    counts = read_tof_maps(build_tof(cspec, neutrons=200_000, seed=1).run(), bins)[MONITOR]
    comparison = compare.mcstas_vs_tof(MONITOR, counts, counts, 1, 200_000, 1.0)
    time_on_source = sc.scalar(cspec.params.t0_offset, unit="s")
    flights = subpulse_flights(comparison, cspec.positions(), time_on_source, "second")
    assert flights.sizes["subpulse"] >= 13
    assert np.all(np.diff(flights.coords["time_at_mc_monitor"].values) > 0)
    energy = flights.coords["from_single_monitor_incident_energy"].to(unit="meV").values
    reference = flights.coords["reference_energy"].to(unit="meV").values
    np.testing.assert_allclose(energy, reference, rtol=0.02)
