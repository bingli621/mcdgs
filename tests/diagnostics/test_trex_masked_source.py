import pytest

from mcdgs.diagnostics.trex.masked_source import (
    LAMBDA_RESOLUTION,
    TIME_RESOLUTION,
    TRexToSampleMonitor,
    define_instrument,
)
from mcdgs.instruments.trex.model import TRex


def test_the_instrument_ends_at_the_sample_monitor_with_finer_bins():
    plain = define_instrument()
    names = list(plain.components)
    assert names == [
        name for name in TRex().components if name not in ("beamstop", "banana_detector")
    ]
    assert names[-1] == "sample_monitor"
    (time_start, time_stop, time_bins), (lambda_start, lambda_stop, lambda_bins) = (
        plain.sample_monitor.bin_limits.values()
    )
    assert (time_stop - time_start) / time_bins == pytest.approx(TIME_RESOLUTION)
    assert (lambda_stop - lambda_start) / lambda_bins == pytest.approx(LAMBDA_RESOLUTION)


def test_each_source_has_a_build_of_its_own():
    plain, masked = define_instrument(), define_instrument(source_masked=True)
    assert isinstance(masked, TRexToSampleMonitor)
    assert (plain.name, masked.name) == ("trex_to_sample_monitor", "trex_to_sample_monitor_masked")
    assert type(plain.source).__name__ == "ESS_butterfly"
    assert type(masked.source).__name__ == "Polygon_ESS_butterfly"
