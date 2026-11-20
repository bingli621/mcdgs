"""McStasScript's diagnostics of T-Rex: the beam along the guides, and the intensity before each component.

:func:`main` runs T-Rex in the white beam, 0.1 to 10 Å, with every chopper parked at 0 Hz on
the slit of the mode, ``rrm`` 0, with :class:`.McStasScriptDiagnostics`, twice, and writes:

- ``trex_guides_diagnostics.png``: the beam after the source, before the first bandwidth chopper,
  at the end of the curved guide, at the end of the straight guide and at the sample position;
  its profile, wavelength against time, and its spectra in time and in wavelength;
- ``trex_guides_transmission.png``: the transmission by wavelength of the sections of guide
  between the points, and from the source and from ``bw1`` to the sample position, from the rays
  of the first run, with :func:`transmission` and :func:`plot_transmission`;
- ``trex_intensity.png``: McStasScript's diagram of the components, with the intensity and the
  number of rays just before each.

They are the figures of the Guide transmission page of T-Rex, and its code that of the how-to
page on McStasScript's diagnostics. Run the module for them::

    pixi run python -m mcdgs.diagnostics.trex.guides_diagnostics --out docs/instruments/trex/figures/guides
"""

from collections.abc import Mapping
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from numpy.typing import NDArray

from mcdgs.diagnostics import cli
from mcdgs.instruments.trex.model import Params, TRex
from mcdgs.mcstas import McStasScriptDiagnostics

NAME = "trex_guides_diagnostics"

# the wavelength bins of the transmission: 0.1 Å wide, over the white beam, 0.1 to 10 Å
WAVELENGTH_BINS = np.linspace(0.1, 10.0, 100)
# a bin with fewer rays at the end of a section has no transmission: too few to be read
MIN_RAYS = 20
# the points by the name the figure gives them
POINT_NAMES = {
    "after source": "source",
    "before bw1": "bw1",
    "after guide_29": "end of the curve",
    "before sample_pos": "sample position",
}
# the sections of guide between consecutive points, and the whole way from two of them; the end
# of the straight guide, after guide_50a, is left out, as the beam reaches the sample unchanged
SECTIONS = (
    ("after source", "before bw1"),
    ("before bw1", "after guide_29"),
    ("after guide_29", "before sample_pos"),
)
OVERALL = (
    ("after source", "before sample_pos"),
    ("before bw1", "before sample_pos"),
)


# --8<-- [start:instrument]
def define_instrument() -> TRex:
    """Return the instrument of this diagnostic: T-Rex in the white beam, its choppers parked."""
    # the white beam, 0.1 to 10 Å, every chopper parked on its wide slit (FF); the rest as default
    params = Params(lambda0=None, rrm=0)
    return TRex(params=params)


# --8<-- [end:instrument]


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

    Returns
    -------
    :
        Two panels on the same wavelength axis: the transmission of each section of
        ``SECTIONS``, on a linear scale, and of each way of ``OVERALL``, on a log scale.
    """
    fig, (sections, overall) = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    # one color per way, in the order of the default cycle, across both panels
    colors = iter(plt.rcParams["axes.prop_cycle"].by_key()["color"])
    for ax, pairs, scale in ((sections, SECTIONS, "linear"), (overall, OVERALL, "log")):
        for start, end in pairs:
            ax.stairs(
                transmission(intensity, rays, start, end),
                bins,
                label=f"{POINT_NAMES[start]} → {POINT_NAMES[end]}",
                color=next(colors),
                linewidth=1.5,
            )
        ax.set(xlabel="Wavelength [Å]", ylabel="Transmission", yscale=scale)
        ax.set_xlim(bins[0], bins[-1])
        ax.grid(alpha=0.3)
        ax.legend()
    sections.set_ylim(0, 1)
    sections.set_title("Sections of guide")
    overall.set_title("To the sample position")
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def main(
    out: Path | None = None,
    ncount: float = 1e6,
    analysis_ncount: float = 1e7,
    seed: int = 1,
) -> list[Path]:
    """Run T-Rex with McStasScript's diagnostics, then for its diagram, and save both figures.

    The instrument is :class:`.TRex` in the white beam with its choppers parked, ``rrm`` 0, so
    the beam reaches the sample position unchopped, but for what the parked disks block.
    The runs go to the run folders ``trex_guides_diagnostics`` and
    ``trex_guides_diagnostics_analysis``, replaced on every run.

    Parameters
    ----------
    out:
        Folder of the figures. If ``None``, ``trex_guides_diagnostics`` in the results folder of
        :meth:`.Workspace.find`.
    ncount:
        Number of McStas neutrons of the run with the diagnostics, which keeps every ray at
        each point.
    analysis_ncount:
        Number of McStas neutrons of the run of the diagram.
    seed:
        Seed of the random numbers.

    Returns
    -------
    :
        The paths of the diagnostics, of the transmission and of the diagram.
    """
    folder = cli.output_folder(out, NAME)

    # --8<-- [start:diagnostics]
    diag = McStasScriptDiagnostics(define_instrument())

    diag.add_point(after="source")
    diag.add_point(before="bw1")  # the first bandwidth chopper
    diag.add_point(after="guide_29")  # the end of the curved guide
    diag.add_point(after="guide_50a")  # the end of the straight guide
    diag.add_point(before="sample_pos")  # the sample position

    diag.add_view("x", "y", bins=100)  # the beam profile
    diag.add_view("t", "l", bins=100)  # wavelength against time
    diag.add_view("t")  # the time spectrum
    diag.add_view("l")  # the wavelength spectrum

    run = diag.run(NAME, ncount=ncount, seed=seed)
    # --8<-- [end:diagnostics]
    print(f"Diagnostics in {run.folder}")

    # --8<-- [start:plot]
    fig = diag.plot()
    diagnostics = cli.save_figure(fig, folder / f"{NAME}.png", dpi=100)
    # --8<-- [end:plot]

    # --8<-- [start:transmission]
    intensity = diag.histogram("l", WAVELENGTH_BINS)  # n/s per bin, at each point
    rays = diag.histogram("l", WAVELENGTH_BINS, weighted=False)
    title = diag.simulation.instrument.running_conditions()
    fig = plot_transmission(intensity, rays, WAVELENGTH_BINS, title)
    transmission_path = cli.save_figure(fig, folder / "trex_guides_transmission.png", dpi=100)
    # --8<-- [end:transmission]

    # --8<-- [start:analysis]
    run, fig = diag.run_analysis(f"{NAME}_analysis", ncount=analysis_ncount, seed=seed)
    diagram = cli.save_figure(fig, folder / "trex_intensity.png", dpi=100)
    # --8<-- [end:analysis]
    return [diagnostics, transmission_path, diagram]


if __name__ == "__main__":
    main(out=cli.parse_out(__doc__))
