from mcdgs.diagnostics.trex import elastic_resolution_incoherent_long as long
from mcdgs.diagnostics.trex import elastic_resolution_incoherent_middle as middle
from mcdgs.diagnostics.trex import elastic_resolution_incoherent_short as short


def test_the_middle_instrument_is_the_short_one_at_2_angstrom_with_the_bender_out():
    instrument = middle.define_instrument("RR")
    assert instrument.params.lambda0 == 2.0
    assert not instrument.params.bender_in
    assert instrument.params.mode == "RR"
    assert list(instrument.components) == list(short.define_instrument("RR").components)


def test_the_middle_runs_do_not_overwrite_the_short_or_long_ones():
    for mode in short.MODES:
        assert set(middle.run_names(mode)).isdisjoint(short.run_names(mode))
        assert set(middle.run_names(mode)).isdisjoint(long.run_names(mode))
