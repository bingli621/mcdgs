from mcdgs.diagnostics.trex.overview import along_beam_range
from mcdgs.instruments.trex.model import TRex


def test_the_view_along_the_beam_leaves_out_the_source_and_the_banana_detector():
    names = list(TRex().components)
    first, stop = along_beam_range(TRex())
    assert names[first - 1] == "source"
    assert names[stop] == "banana_detector"
    assert names[first:stop][0] == "source_monitor"
    assert names[first:stop][-1] == "beamstop"
