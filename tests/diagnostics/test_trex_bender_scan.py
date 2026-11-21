import numpy as np

from mcdgs.diagnostics.trex.bender_scan import (
    BAND_LIMITS,
    BAND_WIDTH,
    define_instrument,
)
from mcdgs.reduction.binning import monitor_bins


def test_every_edge_of_a_band_is_an_edge_of_a_bin_of_the_beam_monitor():
    edges = monitor_bins(define_instrument().beam_monitor.bin_limits)["wavelength"]
    assert define_instrument().beam_monitor.bin_limits["lambda"][:2] == BAND_LIMITS
    for edge in np.arange(BAND_LIMITS[0], BAND_LIMITS[1] + BAND_WIDTH / 2, BAND_WIDTH):
        assert np.isclose(edges.values, edge, rtol=0, atol=1e-12).any()


def test_the_instrument_of_the_scan_is_the_extraction_with_the_bender_in():
    from mcdgs.diagnostics.trex.bender_scan import TRexBenderScan
    from mcdgs.diagnostics.trex.extraction import Extraction, ExtractionParams

    instrument = define_instrument()
    assert isinstance(instrument, TRexBenderScan) and instrument.name == "trex_bender_scan"
    params = ExtractionParams(bender_in=True)
    extraction = Extraction(params)
    assert list(instrument.components) == list(extraction.components)
    # all but the beam monitor, which covers only the bands
    others = {name: c for name, c in instrument.components.items() if name != "beam_monitor"}
    assert others == {n: c for n, c in extraction.components.items() if n != "beam_monitor"}
    assert instrument.instr_params == extraction.instr_params  # bender_in = 1
