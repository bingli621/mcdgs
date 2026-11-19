"""McStasScript's diagnostics of T-Rex: the beam along the guides, and the intensity before each component.

:func:`main` runs T-Rex at its default control parameters with :class:`.McStasScriptDiagnostics`,
twice, and writes:

- ``trex_beam_diagnostics.png``: the beam after the source, before the first bandwidth chopper
  and at the end of the curved guide; its profile, wavelength against x, and its spectrum;
- ``trex_intensity.png``: McStasScript's diagram of the components, with the intensity and the
  number of rays just before each.

They are the figures of the how-to page on McStasScript's diagnostics. Run the module for them::

    pixi run python -m mcdgs.diagnostics.trex.beam_diagnostics --out docs/how-to/figures/mcstasscript-diagnostics
"""

from pathlib import Path

from mcdgs.diagnostics import cli
from mcdgs.instruments.trex.model import Params, TRex
from mcdgs.mcstas import McStasScriptDiagnostics

NAME = "trex_beam_diagnostics"


# --8<-- [start:instrument]
def define_instrument() -> TRex:
    """Return the instrument of this diagnostic: :class:`.TRex` at its default parameters."""
    params = Params()  # the default control parameters
    return TRex(params=params)


# --8<-- [end:instrument]


def main(
    out: Path | None = None,
    ncount: float = 1e6,
    analysis_ncount: float = 1e7,
    seed: int = 1,
) -> list[Path]:
    """Run T-Rex with McStasScript's diagnostics, then for its diagram, and save both figures.

    The runs go to the run folders ``trex_beam_diagnostics`` and
    ``trex_beam_diagnostics_analysis``, replaced on every run.

    Parameters
    ----------
    out:
        Folder of the figures. If ``None``, ``trex_beam_diagnostics`` in the results folder of
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
        The paths of the diagnostics and of the diagram.
    """
    folder = cli.output_folder(out, NAME)

    # --8<-- [start:diagnostics]
    diag = McStasScriptDiagnostics(define_instrument())

    diag.add_point(after="source")
    diag.add_point(before="bw1")  # the first bandwidth chopper
    diag.add_point(after="guide_29")  # the end of the curved guide

    diag.add_view("x", "y", bins=100)  # the beam profile
    diag.add_view("t", "l", bins=100)  # wavelength against x
    # diag.add_view("l")  # the spectrum

    run = diag.run(NAME, ncount=ncount, seed=seed)
    # --8<-- [end:diagnostics]
    print(f"Diagnostics in {run.folder}")

    # --8<-- [start:plot]
    fig = diag.plot()
    diagnostics = cli.save_figure(fig, folder / f"{NAME}.png", dpi=100)
    # --8<-- [end:plot]

    # --8<-- [start:analysis]
    run, fig = diag.run_analysis(f"{NAME}_analysis", ncount=analysis_ncount, seed=seed)
    diagram = cli.save_figure(fig, folder / "trex_intensity.png", dpi=100)
    # --8<-- [end:analysis]
    return [diagnostics, diagram]


if __name__ == "__main__":
    main(out=cli.parse_out(__doc__))
