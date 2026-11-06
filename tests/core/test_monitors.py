import pytest

from mcdgs.core.components.monitors import monitor_time_limits


def test_the_time_window_runs_from_the_fastest_to_the_slowest_leaving_at_the_end():
    # 1 Å flies 3956.034 m/s and 2 Å half that, over 3956.034 m: 1 s and 2 s, plus 0.5 s
    assert monitor_time_limits(3956.034006119441, (1.0, 2.0), 0.5) == pytest.approx((1.0, 2.5))
