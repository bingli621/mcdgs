from mcdgs.diagnostics.trex.extraction import Extraction, ExtractionParams
from mcdgs.instruments.trex.guides.bispectral import BENDER_ROTATION
from mcdgs.instruments.trex.model import TRex

UP_TO_GUIDE_4 = [
    "origin", "source", "arm_nboa", "guide_1", "guide_2", "guide_3",
    "arm_guide_4", "guide_4a", "bender", "guide_4b",
]  # fmt: skip


def test_the_extraction_is_trex_up_to_guide_4():
    extraction = Extraction()
    assert list(extraction.components) == UP_TO_GUIDE_4
    trex = TRex()
    for name, component in extraction.components.items():
        assert component == trex.components[name], name


def test_the_instrument_parameters_follow_the_control_parameters():
    extraction = Extraction()
    assert [(p.name, p.value) for p in extraction.instr_params] == [
        ("bender_in", 0),
        ("b_rot", BENDER_ROTATION),
    ]
    extraction.update_params(bender_in=True, b_rot=-0.5)
    assert extraction.params == ExtractionParams(bender_in=True, b_rot=-0.5)
    assert [(p.name, p.value) for p in extraction.instr_params] == [
        ("bender_in", 1),
        ("b_rot", -0.5),
    ]
