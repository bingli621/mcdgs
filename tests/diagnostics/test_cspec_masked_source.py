import pytest

from mcdgs.diagnostics.cspec.masked_source import (
    LAMBDA_RESOLUTION,
    TIME_RESOLUTION,
    CSPECToSampleMonitor,
    define_instrument,
)
from mcdgs.instruments.cspec.model import CSPEC


def test_the_instrument_ends_at_the_sample_monitor_with_finer_bins():
    plain = define_instrument()
    names = list(plain.components)
    assert names == list(CSPEC().components)
    assert names[-1] == "sample_monitor"
    (time_start, time_stop, time_bins), (lambda_start, lambda_stop, lambda_bins) = (
        plain.sample_monitor.bin_limits.values()
    )
    assert (time_stop - time_start) / time_bins == pytest.approx(TIME_RESOLUTION)
    assert (lambda_stop - lambda_start) / lambda_bins == pytest.approx(LAMBDA_RESOLUTION)


def test_each_source_has_a_build_of_its_own():
    plain, masked = define_instrument(), define_instrument(source_masked=True)
    assert isinstance(masked, CSPECToSampleMonitor)
    assert (plain.name, masked.name) == (
        "cspec_to_sample_monitor",
        "cspec_to_sample_monitor_masked",
    )
    assert type(plain.source).__name__ == "ESS_butterfly"
    assert type(masked.source).__name__ == "Polygon_ESS_butterfly"
