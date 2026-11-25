import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from mcdgs.diagnostics.guides_diagnostics import WAVELENGTH_BINS, plot_transmission, transmission
from mcdgs.diagnostics.trex import guides_diagnostics as trex_guides


def fake_histogram(self, variable, bins, weighted=True):
    """Half the beam lost at each point of T-Rex, in every bin, with enough rays everywhere."""
    n = len(bins) - 1
    names = ["after source", "before bw1", "after guide_29", "after guide_50a", "before sample_pos"]
    return {name: np.full(n, (1e10 if weighted else 1e3) * 0.5**i) for i, name in enumerate(names)}


def test_the_transmission_has_no_value_where_too_few_rays_arrive():
    intensity = {"a": np.array([10.0, 10.0, 0.0]), "b": np.array([4.0, 1.0, 1.0])}
    rays = {"a": np.array([100.0, 100.0, 0.0]), "b": np.array([40.0, 5.0, 40.0])}
    result = transmission(intensity, rays, "a", "b")
    assert result[0] == 0.4
    assert np.isnan(result[1:]).all()  # too few rays; nothing at the start


def test_the_transmission_figure_has_the_sections_and_the_whole_way():
    bins = WAVELENGTH_BINS
    intensity = fake_histogram(None, "l", bins)
    rays = fake_histogram(None, "l", bins, weighted=False)
    fig = plot_transmission(
        intensity,
        rays,
        bins,
        "title",
        point_names=trex_guides.POINT_NAMES,
        sections=trex_guides.SECTIONS,
        overall=trex_guides.OVERALL,
    )
    sections, overall = fig.axes
    assert [t.get_text() for t in sections.get_legend().get_texts()] == [
        "source → bw1",
        "bw1 → end of the curve",
        "end of the curve → sample position",
    ]
    assert [t.get_text() for t in overall.get_legend().get_texts()] == [
        "source → sample position",
        "bw1 → sample position",
    ]
    assert overall.get_yscale() == "log" and fig.get_suptitle() == "title"
    plt.close(fig)
