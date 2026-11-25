"""The transmission of the guides by wavelength, from McStasScript's diagnostics, for any instrument.

:func:`transmission` divides the wavelength spectra of two points of
:meth:`.McStasScriptDiagnostics.histogram`, and :func:`plot_transmission` draws it for the
sections of guide between the points, and for the whole way to the sample position. The guides
diagnostic of each instrument gives its own points.
"""

from collections.abc import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from numpy.typing import NDArray

# the wavelength bins of the transmission: 0.1 Å wide, over the white beam, 0.1 to 10 Å
WAVELENGTH_BINS = np.linspace(0.1, 10.0, 100)
# a bin with fewer rays at the end of a section has no transmission: too few to be read
MIN_RAYS = 20


def transmission(
    intensity: Mapping[str, NDArray[np.float64]],
    rays: Mapping[str, NDArray[np.float64]],
    start: str,
    end: str,
) -> NDArray[np.float64]:
    """Compute the transmission from one point to another, in each bin.

    Parameters
    ----------
    intensity:
        The intensity in each bin by point, in n/s, as
        :meth:`.McStasScriptDiagnostics.histogram` gives it.
    rays:
        The number of rays in each bin by point, on the same bins.
    start, end:
        The points, such as ``"before bw1"``.

    Returns
    -------
    :
        The intensity at ``end`` over that at ``start``; NaN where fewer than ``MIN_RAYS``
        rays reach ``end``.

    Examples
    --------
    >>> intensity = {"a": np.array([10.0, 10.0]), "b": np.array([5.0, 1.0])}
    >>> rays = {"a": np.array([100.0, 100.0]), "b": np.array([50.0, 10.0])}
    >>> transmission(intensity, rays, "a", "b")
    array([0.5, nan])
    """
    ratio = np.divide(
        intensity[end],
        intensity[start],
        out=np.full_like(intensity[end], np.nan),
        where=intensity[start] > 0,
    )
    return np.where(rays[end] >= MIN_RAYS, ratio, np.nan)


def plot_transmission(
    intensity: Mapping[str, NDArray[np.float64]],
    rays: Mapping[str, NDArray[np.float64]],
    bins: NDArray[np.float64],
    title: str,
    *,
    point_names: Mapping[str, str],
    sections: Sequence[tuple[str, str]],
    overall: Sequence[tuple[str, str]],
) -> Figure:
    """Draw the transmission by wavelength of the sections of guide, and of the whole way.

    Parameters
    ----------
    intensity, rays:
        The intensity, in n/s, and the number of rays in each wavelength bin by point, as for
        :func:`transmission`.
    bins:
        The edges of the wavelength bins, in Å.
    title:
        The title of the figure, such as the running conditions of the instrument.
    point_names:
        The name the figure gives each point, such as ``"bw1"`` for ``"before bw1"``.
    sections:
        The sections of guide between consecutive points, as (start, end) points.
    overall:
        The ways to the sample position, likewise.

    Returns
    -------
    :
        Two panels on the same wavelength axis: the transmission of each of ``sections``, on a
        linear scale, and of each of ``overall``, on a log scale.
    """
    fig, (sections_ax, overall_ax) = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    # one color per way, in the order of the default cycle, across both panels
    colors = iter(plt.rcParams["axes.prop_cycle"].by_key()["color"])
    for ax, pairs, scale in ((sections_ax, sections, "linear"), (overall_ax, overall, "log")):
        for start, end in pairs:
            ax.stairs(
                transmission(intensity, rays, start, end),
                bins,
                label=f"{point_names[start]} → {point_names[end]}",
                color=next(colors),
                linewidth=1.5,
            )
        ax.set(xlabel="Wavelength [Å]", ylabel="Transmission", yscale=scale)
        ax.set_xlim(bins[0], bins[-1])
        ax.grid(alpha=0.3)
        ax.legend()
    sections_ax.set_ylim(0, 1)
    sections_ax.set_title("Sections of guide")
    overall_ax.set_title("To the sample position")
    fig.suptitle(title)
    fig.tight_layout()
    return fig
