import math
from dataclasses import dataclass

import yaml

from mcdgs.core.base import INSTRUMENTS, Instrument
from mcdgs.core.components import Progress_bar, Slit
from mcdgs.core.params import InstrParam, param
from mcdgs.instruments import TRex
from mcdgs.instruments.trex import PARAMS_FILE
from mcdgs.io import load, read_yaml, save, to_dict, write_yaml
from mcdgs.samples import la2ni7


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


def test_params_yaml_is_the_default_trex_saved():
    # if this fails after a change to the model: save(TRex(), PARAMS_FILE) and review the diff
    assert _close(read_yaml(PARAMS_FILE), to_dict(TRex()))
    assert _close(to_dict(load(PARAMS_FILE)), to_dict(TRex()))


def test_write_then_read_gives_the_same_dict(tmp_path):
    data = {"name": "x", "a": {"b": [1, 2.5, None], "unit": "Å"}}
    write_yaml(data, tmp_path / "x.yaml")
    assert read_yaml(tmp_path / "x.yaml") == data
    assert "unit: Å" in (tmp_path / "x.yaml").read_text(encoding="utf-8")  # not escaped


def test_an_empty_file_is_an_empty_dict(tmp_path):
    (tmp_path / "empty.yaml").touch()
    assert read_yaml(tmp_path / "empty.yaml") == {}


def test_a_number_in_exponent_form_is_a_number_with_or_without_a_dot(tmp_path):
    (tmp_path / "n.yaml").write_text("a: 1e-3\nb: 2E3\nc: -1.5e+2\nd: 1e\ne: 0x10\n")
    assert read_yaml(tmp_path / "n.yaml") == {
        "a": 0.001,
        "b": 2000.0,
        "c": -150.0,
        "d": "1e",
        "e": 16,
    }
    assert yaml.safe_load("a: 1e-3") == {"a": "1e-3"}  # PyYAML's own loader is not changed


def test_save_then_load_gives_the_instrument_back(tmp_path):
    trex = TRex(sample=la2ni7(mosaic=60))
    trex.update_params(lambda0=4.0, bender_in=True, mode="RR")
    save(trex, tmp_path / "trex.yaml")
    assert load(tmp_path / "trex.yaml") == trex  # the sample comes from the file


def test_an_entry_is_written_on_one_line_however_long(tmp_path):
    save(TRex(), tmp_path / "trex.yaml")
    text = (tmp_path / "trex.yaml").read_text(encoding="utf-8")
    assert "  t0_offset: {value: 0.0017, unit: s, comment: time added to choppers" in text
    assert "  at: [0, 0, 0]\n" in text
    assert not any(line.startswith("      ") for line in text.splitlines())  # never wrapped


def test_a_file_is_laid_over_the_code(tmp_path):
    write_yaml(
        {"instrument": "TRex", "params": {"lambda0": 4.0}, "slit_mc": {"xwidth": 0.02}},
        tmp_path / "over.yaml",
    )
    trex = load(tmp_path / "over.yaml")
    assert (trex.params.lambda0, trex.slit_mc.xwidth) == (4.0, 0.02)
    assert trex.slit_mc.yheight == TRex().slit_mc.yheight  # left out: as made


def test_a_sample_a_file_adds_is_built_where_the_file_puts_it(tmp_path):
    (tmp_path / "sample.yaml").write_text(
        "instrument: TRex\n"
        "sample_pos: {}\n"
        "sample:\n"
        "  component: Single_crystal\n"
        "  reflections: Al.lau\n"
        "  lattice_params: [4.0495, 4.0495, 4.0495, 90, 90, 90]\n"
        "  along_beam: [0, 0, 1]\n"
        "  in_plane: [1, 0, 0]\n"
        "  relative: sample_pos\n"
    )
    trex = load(tmp_path / "sample.yaml")
    assert list(trex.components)[-2:] == ["sample_pos", "sample"]
    assert trex.sample.relative == "sample_pos"


@dataclass(frozen=True)
class _TinyParams:
    gap: float = param(0.01, "m", "slit opening")


class _Tiny(Instrument):
    """Another instrument: io/ works with it unchanged."""

    params_class = _TinyParams
    name = "tiny"

    def __init__(self, params=None):
        super().__init__(params)
        self.origin = Progress_bar(at=[0, 0, 0])
        self.slit = Slit(xwidth="gap", yheight=0.02, at=[0, 0, 1.0], relative="origin")

    @property
    def instr_params(self):
        return [InstrParam("gap", self.params.gap, "m", "slit opening")]


def test_another_instrument_goes_through_the_same_io(tmp_path):
    assert INSTRUMENTS["_Tiny"] is _Tiny
    tiny = _Tiny(_TinyParams(gap=0.05))
    save(tiny, tmp_path / "tiny.yaml")
    assert read_yaml(tmp_path / "tiny.yaml")["instrument"] == "_Tiny"
    assert load(tmp_path / "tiny.yaml") == tiny  # the file names its class
    write_yaml({"instrument": "_Tiny", "slit": {"yheight": 0.03}}, tmp_path / "over.yaml")
    assert load(tmp_path / "over.yaml").slit.yheight == 0.03  # a partial file names it too
