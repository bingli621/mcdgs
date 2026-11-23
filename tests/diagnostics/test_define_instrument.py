import importlib
import inspect
import pkgutil

import pytest

import mcdgs.diagnostics.cspec as cspec_diagnostics
import mcdgs.diagnostics.trex as trex_diagnostics
from mcdgs.core import Instrument

SHARED = {"mcdgs.diagnostics.trex.sample_runs"}  # shared by the diagnostics, not one itself
MODULES = sorted(
    name
    for package in (trex_diagnostics, cspec_diagnostics)
    for info in pkgutil.iter_modules(package.__path__)
    if not info.ispkg and (name := f"{package.__name__}.{info.name}") not in SHARED
)


def test_every_diagnostic_is_listed():
    trex = [m for m in MODULES if m.startswith(trex_diagnostics.__name__)]
    cspec = [m for m in MODULES if m.startswith(cspec_diagnostics.__name__)]
    assert (len(trex), len(cspec)) == (16, 2)


@pytest.mark.parametrize("name", MODULES)
def test_each_diagnostic_defines_its_instrument_in_one_function(name):
    module = importlib.import_module(name)
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
