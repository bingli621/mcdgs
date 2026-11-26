"""The S(Q, E) of a powder magnon in every sub-pulse, drawn for any instrument.

The powder-magnon diagnostic of each instrument writes the S(|Q|, E) of its powder with
:func:`.write_magnon_sqw`, runs McStas, and reduces its detector to the S(Q, E) of each sub-pulse
with :func:`.q_e_maps`, as for vanadium. Each chooses a magnon in reach of its sub-pulses. This
module does the rest:

- :func:`crop` keeps the bins of the maps within a range of E, around the band.
- :func:`plot_magnon_maps` draws the maps of :func:`.plot_q_e_maps` with the edges of the band of
  the model, :func:`.magnon_band`, over them: how far the intensity spreads past them is the
  resolution of each sub-pulse.
"""

import numpy as np
import scipp as sc
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from mcdgs.diagnostics.elastic_resolution import plot_q_e_maps
from mcdgs.samples.inelastic.powder import MAGNON_A, magnon_band

EDGE_COLOR = "white"  # of the edges of the band of the model, over the maps


def crop(maps: sc.DataArray, energy_range: tuple[float, float]) -> sc.DataArray:
    """Keep the bins of E within ``energy_range``, in meV, of the maps of :func:`.q_e_maps`.

    Examples
    --------
    >>> edges = sc.linspace("energy_transfer", -10.0, 10.0, 21, unit="meV")
    >>> maps = sc.DataArray(
    ...     sc.zeros(dims=["energy_transfer"], shape=[20]), coords={"energy_transfer": edges}
    ... )
    >>> cropped = crop(maps, (-2.0, 3.0))
    >>> cropped.coords["energy_transfer"].values[[0, -1]].tolist()
    [-2.0, 3.0]
    """
    low, high = (sc.scalar(value, unit="meV") for value in energy_range)
    unit = maps.coords["energy_transfer"].unit
    return maps["energy_transfer", low.to(unit=unit) : high.to(unit=unit)]


def plot_magnon_maps(
    maps: sc.DataArray,
    title: str,
    energy_range: tuple[float, float],
    js: float,
    gap: float,
    a: float = MAGNON_A,
) -> Figure:
    """Draw the S(Q, E) of every sub-pulse, with the edges of the band of the model over it.

    The maps are those of :func:`.plot_q_e_maps`, within ``energy_range``: one panel per
    sub-pulse, from the highest Ei, on one color scale, white where the sub-pulse cannot reach.
    Over each, the lowest and the highest energy of the magnon at each |Q|, :func:`.magnon_band`,
    are dashed lines: the band of the powder has sharp edges, and how far the intensity spreads
    past them is the resolution of the sub-pulse.

    Parameters
    ----------
    maps:
        The maps of :func:`.q_e_maps`.
    title:
        The title of the figure.
    energy_range:
        The E of the maps, in meV, as :func:`crop`.
    js, gap, a:
        The magnon of the powder, as :func:`.magnon_energies`.

    Returns
    -------
    :
        The figure.
    """
    maps = crop(maps, energy_range)
    fig = plot_q_e_maps(maps, title)
    q_edges = maps.coords["Q"].to(unit="1/Å").values
    q = np.linspace(q_edges[0], q_edges[-1], 400)
    low, high = magnon_band(q, a=a, js=js, gap=gap)
    panels = [ax for ax in fig.axes[: maps.sizes["subpulse"]] if ax.get_visible()]
    for ax in panels:
        for edge in (low, high):
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
