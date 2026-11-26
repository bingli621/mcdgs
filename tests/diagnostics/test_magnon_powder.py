import numpy as np
import scipp as sc

from mcdgs.diagnostics import magnon_powder as magnon


def test_the_figure_has_a_panel_per_sub_pulse_with_the_band_over_it():
    n = 15
    edges = sc.linspace("energy_transfer", -10.0, 30.0, 401, unit="meV")
    q = sc.linspace("Q", 0.0, 6.0, 61, unit="1/Å")
    maps = sc.DataArray(
        sc.array(
            dims=["subpulse", "energy_transfer", "Q"],
            values=np.random.default_rng(1).random((n, 400, 60)),
        ),
        coords={
            "energy_transfer": edges,
            "Q": q,
            "incident_energy": sc.linspace("subpulse", 47.0, 11.0, n, unit="meV"),
        },
    )
    fig = magnon.plot_magnon_maps(maps, "test", js=15 / 24, gap=5.0)
    panels = [ax for ax in fig.axes if ax.get_title().split(":")[0].isdigit()]
    assert len(panels) == n
    assert all(len(ax.lines) == 4 for ax in panels)  # the two edges of the band, on both sides
