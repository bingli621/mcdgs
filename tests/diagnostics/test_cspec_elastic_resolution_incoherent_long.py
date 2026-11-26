from mcdgs.diagnostics.cspec import elastic_resolution_incoherent_long as long
from mcdgs.diagnostics.cspec import elastic_resolution_incoherent_short as short


def test_the_long_instrument_is_the_short_one_at_8_angstrom():
    instrument = long.define_instrument()
    assert instrument.params.lambda0 == 8.0
    assert list(instrument.components) == list(short.define_instrument().components)


def test_the_long_runs_do_not_overwrite_the_short_ones():
    assert set(long.run_names()).isdisjoint(short.run_names())
