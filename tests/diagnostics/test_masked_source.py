from dataclasses import dataclass

import pytest

from mcdgs.core import INSTRUMENTS, Instrument
from mcdgs.core.components import Arm
from mcdgs.diagnostics.masked_source import run_source


@dataclass(frozen=True)
class NoParams:
    """No control parameters."""


def test_run_source_needs_a_time_and_wavelength_monitor():
    class NoMonitor(Instrument):
        """A source and, where the sample monitor should be, an arm."""

        params_class = NoParams
        name = "no_monitor"

        def __init__(self) -> None:
            super().__init__()
            self.source = Arm()
            self.sample_monitor = Arm()

    try:
        with pytest.raises(TypeError, match="sample_monitor is a Arm, not a Monitor_nD"):
            run_source(NoMonitor(), "never_run", ncount=1, seed=1)
    finally:
        del INSTRUMENTS["NoMonitor"]
