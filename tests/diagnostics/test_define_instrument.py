import importlib
import inspect
import pkgutil

import pytest

import mcdgs.diagnostics.trex as trex_diagnostics
from mcdgs.core import Instrument

MODULES = sorted(
    info.name for info in pkgutil.iter_modules(trex_diagnostics.__path__) if not info.ispkg
)


def test_every_diagnostic_of_t_rex_is_listed():
    assert len(MODULES) == 13


@pytest.mark.parametrize("name", MODULES)
def test_each_diagnostic_defines_its_instrument_in_one_function(name):
    module = importlib.import_module(f"mcdgs.diagnostics.trex.{name}")
    define_instrument = module.define_instrument
    assert define_instrument.__module__ == module.__name__  # its own, not imported
    instrument = define_instrument()  # every argument has a default
    assert isinstance(instrument, Instrument)
    source = inspect.getsource(module)
    assert source.count("# --8<-- [start:instrument]") == 1
    assert source.count("# --8<-- [end:instrument]") == 1
    start = source.index("# --8<-- [start:instrument]")
    assert (
        start < source.index("def define_instrument(") < source.index("# --8<-- [end:instrument]")
    )
