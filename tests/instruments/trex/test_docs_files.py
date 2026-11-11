"""The files of the documentation are committed: check they match the code.

If one fails, run `pixi run docs-diagnostics` and commit what it changes.
"""

from pathlib import Path

from mcdgs.instruments.trex.model import TRex
from mcdgs.mcstas.builder import build_mcstas

DOCS_INSTR = (
    Path(__file__).resolve().parents[3] / "docs" / "instruments" / "trex" / "figures" / "overview"
)


def _without_date(text: str) -> list[str]:
    """The lines of an .instr file, less the date McStasScript writes in its header."""
    return [line for line in text.splitlines() if not line.startswith("* Date:")]


def test_the_instr_file_of_the_docs_matches_the_default_model(tmp_path):
    trex = TRex()
    build_mcstas(trex, input_path=tmp_path).write_full_instrument()
    fresh = (tmp_path / f"{trex.name}.instr").read_text()
    committed = (DOCS_INSTR / f"{trex.name}.instr").read_text()
    assert _without_date(fresh) == _without_date(committed), "run `pixi run docs-diagnostics`"
