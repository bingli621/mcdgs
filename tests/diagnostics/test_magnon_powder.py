import matplotlib

matplotlib.use("Agg")

import numpy as np
import scipp as sc

from mcdgs.diagnostics import powder_magnon as magnon
from mcdgs.diagnostics.cspec import powder_magnon as cspec
from mcdgs.diagnostics.trex import powder_magnon as trex
from mcdgs.samples import write_magnon_sqw
from mcdgs.tof import beamline


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


def test_the_time_distance_diagram_runs_from_the_monitor_before_the_sample_to_the_detector(
    tmp_path,
):
    sqw = write_magnon_sqw(tmp_path / "m.sqw", q_max=2.0, q_step=0.5, directions=1000)
    for diagnostic, instrument in (
        (trex, trex.TRexMagnonPowder(sqw_file=sqw)),
        (cspec, cspec.CSPECMagnonPowder(sqw_file=sqw)),
    ):
        fig = magnon.plot_time_distance(
            instrument, diagnostic.TIME_DISTANCE_MONITOR, "test", neutrons=100_000
        )
        distances = {name: d for name, _, d in beamline(instrument)[1]}
        low, high = fig.axes[0].get_ylim()
        assert fig.axes[0].get_title().startswith("test\n")
        assert low < distances[diagnostic.TIME_DISTANCE_MONITOR] < distances["sample_pos"] < high
        assert distances[diagnostic.TIME_DISTANCE_MONITOR] - low < 0.1 * (high - low)
