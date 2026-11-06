from dataclasses import replace

from mcdgs.diagnostics.trex_01_source import TRexSource
from mcdgs.diagnostics.trex_02_extraction import Extraction, ExtractionParams
from mcdgs.instruments.trex.guides.bispectral import BENDER_ROTATION
from mcdgs.instruments.trex.model import TRex

UP_TO_THE_BEAM_MONITOR = [
    "origin", "source", "source_monitor", "arm_nboa", "guide_1", "guide_2", "guide_3",
    "arm_guide_4", "guide_4a", "bender", "guide_4b", "beam_monitor",
]  # fmt: skip


def test_the_extraction_is_trex_up_to_the_beam_monitor():
    extraction = Extraction()
    assert list(extraction.components) == UP_TO_THE_BEAM_MONITOR
    trex = TRex()
    for name, component in extraction.components.items():
        if name not in ("source", "source_monitor", "beam_monitor"):
            assert component == trex.components[name], name
    # the extraction keeps the widest band, which only the source and the beam monitor hold
    assert extraction.source == replace(trex.source, Lmin=0.1, Lmax=10.0)
    assert extraction.source_monitor == TRexSource().source_monitor
    assert extraction.beam_monitor == replace(
        trex.beam_monitor, options=extraction.beam_monitor.options
    )


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


def test_the_source_monitor_is_the_entrance_of_guide_1():
    extraction = Extraction()
    positions = extraction.positions()
    assert positions["source_monitor"] == positions["guide_1"]
    assert extraction.source_monitor.rotated == extraction.components["arm_nboa"].rotated
    assert (extraction.source_monitor.xwidth, extraction.source_monitor.yheight) == (
        extraction.guide_1.w1,
        extraction.guide_1.h1,
    )
