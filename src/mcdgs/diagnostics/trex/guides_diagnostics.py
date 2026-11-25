"""McStasScript's diagnostics of T-Rex: the beam along the guides, and the intensity before each component.

:func:`main` runs T-Rex in the white beam, 0.1 to 10 Å, with every chopper parked at 0 Hz on
the slit of the mode, ``rrm`` 0, with :class:`.McStasScriptDiagnostics`, twice, and writes:

- ``trex_guides_diagnostics.png``: the beam after the source, before the first bandwidth chopper,
  at the end of the curved guide, at the end of the straight guide and at the sample position;
  its profile, wavelength against time, and its spectra in time and in wavelength;
- ``trex_guides_transmission.png``: the transmission by wavelength of the sections of guide
  between the points, and from the source and from ``bw1`` to the sample position, from the rays
  of the first run, with :func:`.transmission` and :func:`.plot_transmission`;
- ``trex_intensity.png``: McStasScript's diagram of the components, with the intensity and the
  number of rays just before each.

They are the figures of the Guide transmission page of T-Rex, and its code that of the how-to
page on McStasScript's diagnostics. Run the module for them::

    pixi run python -m mcdgs.diagnostics.trex.guides_diagnostics --out docs/instruments/trex/figures/guides
"""

from pathlib import Path

from mcdgs.diagnostics import cli
from mcdgs.diagnostics.guides_diagnostics import WAVELENGTH_BINS, plot_transmission
from mcdgs.instruments.trex.model import Params, TRex
from mcdgs.mcstas import McStasScriptDiagnostics

NAME = "trex_guides_diagnostics"
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
    fig = plot_transmission(
        intensity,
        rays,
        WAVELENGTH_BINS,
        title,
        point_names=POINT_NAMES,
        sections=SECTIONS,
        overall=OVERALL,
    )
    transmission_path = cli.save_figure(fig, folder / "trex_guides_transmission.png", dpi=100)
    # --8<-- [end:transmission]

    # --8<-- [start:analysis]
    run, fig = diag.run_analysis(f"{NAME}_analysis", ncount=analysis_ncount, seed=seed)
    diagram = cli.save_figure(fig, folder / "trex_intensity.png", dpi=100)
    # --8<-- [end:analysis]
    return [diagnostics, transmission_path, diagram]


if __name__ == "__main__":
    main(out=cli.parse_out(__doc__))
