import matplotlib
import pytest

matplotlib.use("Agg")

from matplotlib.figure import Figure

from mcdgs.instruments.trex.model import TRex
from mcdgs.tof.acceptance import chopper_cascade


@pytest.mark.parametrize("mode", ["FF", "RR"])
def test_trex_cascade_accepts_a_band_around_lambda0(mode):
    trex = TRex()
    trex.update_params(mode=mode)
    frames, names = chopper_cascade(trex)
    assert names == ["source", "bw1", "bw2", "ps1", "ps2", "mc1", "mc2"]
    accepted = frames.frames[-1].subframes
    lambda0 = trex.params.lambda0
    assert any(s.wavelength.min().value < lambda0 < s.wavelength.max().value for s in accepted)


def test_trex_cascade_splits_the_band_by_repetition_rate_multiplication():
    trex = TRex()
    frames, _ = chopper_cascade(trex)
    # the pulse-shaping choppers turn at 13.5 times the source frequency: no whole multiple
    centers = sorted(
        {round(float(s.wavelength.mean().value), 2) for s in frames.frames[-1].subframes}
    )
    assert len(centers) > 10
    l_min, l_max = trex.source.Lmin, trex.source.Lmax
    assert all(l_min <= c <= l_max for c in centers)


def test_cascade_takes_a_wider_band_than_the_source():
    frames, _ = chopper_cascade(TRex(), wavelength_min=0.1, wavelength_max=10)
    source = frames.frames[0].subframes[0]
    assert source.wavelength.min().value == pytest.approx(0.1)
    assert source.wavelength.max().value == pytest.approx(10.0)
    # bw1 and bw2 keep the T-Rex band, about 1.81 Å wide, out of 9.9 Å
    after_bw2 = frames.frames[2].subframes
    assert max(s.wavelength.max().value for s in after_bw2 if s.wavelength.min().value < 5) < 3.5


def test_trex_draws_one_legend_entry_per_frame():
    figure = TRex().acceptance_diagram()
    assert isinstance(figure, Figure)
    (ax,) = figure.axes
    labels = [text.get_text() for text in ax.get_legend().get_texts()]
    assert labels[0] == "source"
    assert [label.split()[0] for label in labels[1:]] == ["bw1", "bw2", "ps1", "ps2", "mc1", "mc2"]


def test_the_title_is_the_running_conditions():
    trex = TRex()
    trex.update_params(mode="RR", rrm=12)
    assert trex.acceptance_diagram().axes[0].get_title() == trex.running_conditions()


def test_cascade_takes_the_times_and_pulses_given():
    frames, _ = chopper_cascade(TRex(), time_min=-0.001, time_max=0.005, n_pulses=2)
    pulses = frames.frames[0].subframes
    assert len(pulses) == 2
    first, second = (s.time.to(unit="s").values for s in pulses)
    assert (first.min(), first.max()) == pytest.approx((-0.001, 0.005))
    assert second.min() - first.min() == pytest.approx(1 / 14)
    # the openings reach the slowest neutron of the second pulse
    assert any(s.time.min().value > 1 / 14 for s in frames.frames[-1].subframes)


@pytest.mark.parametrize("kwargs", [{"time_min": 0.002, "time_max": 0.001}, {"n_pulses": 0}])
def test_cascade_refuses_an_empty_pulse(kwargs):
    with pytest.raises(ValueError):
        chopper_cascade(TRex(), **kwargs)
