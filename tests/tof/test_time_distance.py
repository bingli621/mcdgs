import matplotlib

matplotlib.use("Agg")

from matplotlib.figure import Figure

from mcdgs.instruments.trex.model import TRex


def test_the_title_is_the_running_conditions():
    trex = TRex()
    trex.update_params(mode="RR", rrm=12)
    figure = trex.time_distance_diagram(10_000, seed=1)
    assert isinstance(figure, Figure)
    assert figure.axes[0].get_title() == trex.running_conditions()


def test_draws_the_pulses_given():
    figure = TRex().time_distance_diagram(10_000, seed=1, n_pulses=2, blocked_rays=0)
    labels = [text.get_text() for text in figure.axes[0].texts]
    assert labels.count("Pulse") == 2
