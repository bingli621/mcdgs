from mcdgs.diagnostics.trex import elastic_resolution_incoherent_long as long
from mcdgs.diagnostics.trex import elastic_resolution_incoherent_short as short


def test_the_short_instrument_is_at_1_angstrom():
    assert short.define_instrument().params.lambda0 == 1.0


def test_the_long_instrument_is_the_short_one_at_3_angstrom():
    instrument = long.define_instrument("RR")
    assert instrument.params.lambda0 == 3.0
    assert instrument.params.mode == "RR"
    assert list(instrument.components) == list(short.define_instrument("RR").components)


def test_the_long_runs_do_not_overwrite_the_short_ones():
    for mode in short.MODES:
        assert set(long.run_names(mode)).isdisjoint(short.run_names(mode))
