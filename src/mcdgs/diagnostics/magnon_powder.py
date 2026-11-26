"""The S(Q, E) of a powder magnon in every sub-pulse, drawn for any instrument.

The powder-magnon diagnostic of each instrument writes the S(|Q|, E) of its powder with
:func:`.write_magnon_sqw`, runs McStas, and reduces its detector to the S(Q, E) of each sub-pulse
with :func:`.q_e_maps`, as for vanadium. Each chooses a magnon in reach of its sub-pulses. This
module does the rest:

- :func:`plot_magnon_maps` draws the maps of :func:`.plot_q_e_maps` over the whole arrival
  windows, energy gain and energy loss, with the edges of the band of the model,
  :func:`.magnon_band`, over them: how far the intensity spreads past them is the resolution of
  each sub-pulse, and whatever lies elsewhere comes from the instrument or the reduction.
"""

import numpy as np
import scipp as sc
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from mcdgs.diagnostics.elastic_resolution import MAP_PERCENTILE, plot_q_e_maps
from mcdgs.samples.inelastic.powder import MAGNON_A, magnon_band

EDGE_COLOR = "white"  # of the edges of the band of the model, over the maps


def plot_magnon_maps(
    maps: sc.DataArray,
    title: str,
    js: float,
    gap: float,
    a: float = MAGNON_A,
    percentile: float = MAP_PERCENTILE,
) -> Figure:
    """Draw the S(Q, E) of every sub-pulse, with the edges of the band of the model over it.

    The maps are those of :func:`.plot_q_e_maps`, over all their E, from the energy gain to the
    energy loss of the arrival window of the highest Ei: one panel per sub-pulse, from the
    highest Ei, on one color scale, white where the sub-pulse cannot reach.

    Over each, the lowest and the highest energy of the magnon at each |Q|, :func:`.magnon_band`,
    are dashed lines, at positive E, where the powder creates the magnon, and mirrored at negative
    E, where it annihilates it: the band of the powder has sharp edges, and how far the intensity
    spreads past them is the resolution of the sub-pulse. Whatever lies outside both bands comes
    from the instrument or the reduction.

    Parameters
    ----------
    maps:
        The maps of :func:`.q_e_maps`.
    title:
        The title of the figure.
    js, gap, a:
        The magnon of the powder, as :func:`.magnon_energies`.
    percentile:
        The top of the color scale, as :func:`.plot_q_e_maps`.

    Returns
    -------
    :
        The figure.
    """
    fig = plot_q_e_maps(maps, title, percentile=percentile)
    q_edges = maps.coords["Q"].to(unit="1/Å").values
    q = np.linspace(q_edges[0], q_edges[-1], 400)
    low, high = magnon_band(q, a=a, js=js, gap=gap)
    panels = [ax for ax in fig.axes[: maps.sizes["subpulse"]] if ax.get_visible()]
    for ax in panels:
        for edge in (low, high, -low, -high):  # on the loss side, and mirrored on the gain side
            ax.plot(q, edge, color=EDGE_COLOR, lw=0.8, ls="--")
    if panels:
        panels[0].legend(
            [Line2D([], [], color=EDGE_COLOR, lw=0.8, ls="--")],
            ["Edges of the band, without resolution"],
            loc="upper right",
            fontsize="small",
            facecolor="0.3",
            labelcolor=EDGE_COLOR,
        )
    return fig
