"""What every instrument does the same way: its name, its parameters file, its committed
instrument file, its sample and its source."""

import math
from pathlib import Path

import pytest

from mcdgs.core import INSTRUMENTS
from mcdgs.core.components import ESS_butterfly, Polygon_ESS_butterfly
from mcdgs.instruments import CSPEC, TRex, cspec, trex
from mcdgs.io import load, read_yaml, to_dict
from mcdgs.mcstas import build_mcstas
from mcdgs.samples import aluminum, la2ni7, lead

DOCS = Path(__file__).resolve().parents[2] / "docs" / "instruments"

# each instrument with its parameters file and the folder of its pages
CASES = [
    pytest.param(TRex, trex.PARAMS_FILE, "trex", id="TRex"),
    pytest.param(CSPEC, cspec.PARAMS_FILE, "cspec", id="CSPEC"),
]
INSTRUMENT_CLASSES = [pytest.param(TRex, id="TRex"), pytest.param(CSPEC, id="CSPEC")]


def _close(a, b):
    """`a` and `b` are the same nested dict, numbers but bools equal to 1e-12 relative: the math
    library of another platform may round atan or sin in the last bit."""
    if isinstance(a, dict):
        return isinstance(b, dict) and list(a) == list(b) and all(_close(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return isinstance(b, list) and len(a) == len(b) and all(map(_close, a, b))
    numbers = (int, float)
    if isinstance(a, numbers) and isinstance(b, numbers) and bool not in (type(a), type(b)):
        return math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-15)
    return a == b and type(a) is type(b)


def _without_date(text: str) -> list[str]:
    """The lines of an .instr file, less the date McStasScript writes in its header."""
    return [line for line in text.splitlines() if not line.startswith("* Date:")]


@pytest.mark.parametrize("cls", INSTRUMENT_CLASSES)
def test_a_file_can_name_the_instrument(cls):
    assert INSTRUMENTS[cls.name] is cls


@pytest.mark.parametrize(("cls", "params_file", "_"), CASES)
def test_params_yaml_is_the_default_instrument_saved(cls, params_file, _):
    # if this fails after a change to the model: save(cls(), params_file) and review the diff
    assert _close(read_yaml(params_file), to_dict(cls()))
    assert _close(to_dict(load(params_file)), to_dict(cls()))


@pytest.mark.parametrize(("cls", "_", "folder"), CASES)
def test_the_instr_file_of_the_docs_matches_the_default_model(cls, _, folder, tmp_path):
    # the files of the documentation are committed: if this fails, run
    # `pixi run docs-diagnostics` and commit what it changes
    instrument = cls()
    build_mcstas(instrument, input_path=tmp_path).write_full_instrument()
    fresh = (tmp_path / f"{instrument.name}.instr").read_text()
    committed = (DOCS / folder / "figures" / "overview" / f"{instrument.name}.instr").read_text()
    assert _without_date(fresh) == _without_date(committed), "run `pixi run docs-diagnostics`"


@pytest.mark.parametrize("cls", INSTRUMENT_CLASSES)
def test_without_a_sample_there_is_none(cls):
    instrument = cls()
    assert instrument.sample is None
    assert "sample" not in instrument.components


@pytest.mark.parametrize("cls", INSTRUMENT_CLASSES)
def test_a_sample_is_at_the_sample_position(cls):
    instrument = cls(sample=aluminum())
    assert instrument.sample.relative == "sample_pos"
    assert instrument.positions()["sample"] == instrument.positions()["sample_pos"]


@pytest.mark.parametrize("make", [aluminum, la2ni7, lead])
@pytest.mark.parametrize("cls", INSTRUMENT_CLASSES)
def test_every_ready_made_sample_goes_on_the_instrument(cls, make):
    instrument = cls(sample=make())
    assert type(instrument.sample) is type(make())
    assert instrument.sample == make(relative="sample_pos")


@pytest.mark.parametrize("cls", INSTRUMENT_CLASSES)
def test_a_sample_relative_to_the_sample_position_stays_there(cls):
    assert cls(sample=aluminum(relative="sample_pos")).sample.relative == "sample_pos"


@pytest.mark.parametrize("cls", INSTRUMENT_CLASSES)
def test_a_sample_relative_to_another_component_is_refused(cls):
    with pytest.raises(ValueError, match="'source', not to 'sample_pos'"):
        cls(sample=aluminum(relative="source"))


@pytest.mark.parametrize("cls", INSTRUMENT_CLASSES)
def test_the_given_sample_is_not_changed(cls):
    crystal = aluminum()
    instrument = cls(sample=crystal)
    assert instrument.sample is not crystal
    assert (crystal.name, crystal.relative) == (None, None)


@pytest.mark.parametrize("cls", INSTRUMENT_CLASSES)
def test_by_default_the_source_is_not_masked(cls):
    assert type(cls().source) is ESS_butterfly


@pytest.mark.parametrize("masked", [False, True])
@pytest.mark.parametrize("cls", INSTRUMENT_CLASSES)
def test_a_new_lambda0_keeps_the_source_masked_or_not(cls, masked):
    instrument = cls(source_masked=masked)
    instrument.update_params(lambda0=4.0)
    assert isinstance(instrument.source, Polygon_ESS_butterfly) == masked
    assert instrument.source == cls(params=instrument.params, source_masked=masked).source
