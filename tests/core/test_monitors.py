import pytest

from mcdgs.core.components import Monitor_nD
from mcdgs.core.components.monitors import monitor_time_limits


def test_the_time_window_runs_from_the_fastest_to_the_slowest_leaving_at_the_end():
    # 1 Å flies 3956.034 m/s and 2 Å half that, over 3956.034 m: 1 s and 2 s, plus 0.5 s
    assert monitor_time_limits(3956.034006119441, (1.0, 2.0), 0.5) == pytest.approx((1.0, 2.5))


@pytest.mark.parametrize(
    ("monitor", "radius"),
    [
        (Monitor_nD(xwidth=0.1, yheight=0.1, options="lambda limits=[1 2] bins=10"), None),
        (Monitor_nD(radius=3.0, yheight=2.0, options="banana, theta limits=[0 90] bins=9"), 3.0),
        (Monitor_nD(xwidth=6.0, yheight=2.0, options="banana, theta limits=[0 90] bins=9"), 3.0),
        (Monitor_nD(zdepth=4.0, options="sphere"), 2.0),
        (Monitor_nD(radius=1.5, yheight=1.0), 1.5),  # McStas makes it a banana
    ],
)
def test_the_curved_radius_is_that_mcstas_takes(monitor, radius):
    assert monitor.curved_radius == radius
