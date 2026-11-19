"""McStasScript's diagnostics of T-Rex: the beam along the guides, and the intensity before each component.

:func:`main` runs T-Rex in the white beam, 0.1 to 10 Å, with every chopper parked at 0 Hz
(:class:`ParkedTRex`), with :class:`.McStasScriptDiagnostics`, twice, and writes:

- ``trex_guides_diagnostics.png``: the beam after the source, before the first bandwidth chopper,
  at the end of the curved guide, at the end of the straight guide and at the sample position;
  its profile, wavelength against time, and its spectra in time and in wavelength;
- ``trex_intensity.png``: McStasScript's diagram of the components, with the intensity and the
  number of rays just before each.

They are the figures of the how-to page on McStasScript's diagnostics. Run the module for them::

    pixi run python -m mcdgs.diagnostics.trex.guides_diagnostics --out docs/how-to/figures/mcstasscript-diagnostics
"""

from pathlib import Path

from mcdgs.core import InstrParam
from mcdgs.diagnostics import cli
from mcdgs.instruments.trex.choppers import DISKS
from mcdgs.instruments.trex.model import Params, TRex
from mcdgs.mcstas import McStasScriptDiagnostics

NAME = "trex_guides_diagnostics"


# --8<-- [start:instrument]
class ParkedTRex(TRex):
    """:class:`.TRex` with every chopper parked, at 0 Hz.

    ``NXdisk_chopper`` then stands still with its zero mark at ``park_angle``, 0, where each disk
    of T-Rex has an opening on the beam, which passes most of it, unchopped. With no chopper to
    phase, it runs with ``lambda0`` ``None``, the white beam.
    """

    def _chopper_instr_params(self) -> list[InstrParam]:
        """Return the frequency, the delay and the park angle of every chopper, all 0: parked."""
        return [
            p
            for name in DISKS
            for p in (
                InstrParam(f"nu_{name}", 0.0, "Hz", f"rotation frequency of chopper {name}"),
                InstrParam(f"delay_{name}", 0.0, "s", f"delay of chopper {name}"),
                InstrParam(f"park_{name}", 0.0, "deg", f"park angle of chopper {name}"),
            )
        ]


def define_instrument() -> ParkedTRex:
    """Return the instrument of this diagnostic: T-Rex in the white beam, its choppers parked."""
    params = Params(lambda0=None)  # the white beam, 0.1 to 10 Å; the rest as default
    return ParkedTRex(params=params)


# --8<-- [end:instrument]


def main(
    out: Path | None = None,
    ncount: float = 1e6,
    analysis_ncount: float = 1e7,
    seed: int = 1,
) -> list[Path]:
    """Run T-Rex with McStasScript's diagnostics, then for its diagram, and save both figures.

    The instrument is :class:`ParkedTRex` in the white beam, so the beam reaches the sample
    position unchopped, but for what the parked disks block.
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
        The paths of the diagnostics and of the diagram.
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

    # --8<-- [start:analysis]
    run, fig = diag.run_analysis(f"{NAME}_analysis", ncount=analysis_ncount, seed=seed)
    diagram = cli.save_figure(fig, folder / "trex_intensity.png", dpi=100)
    # --8<-- [end:analysis]
    return [diagnostics, diagram]


if __name__ == "__main__":
    main(out=cli.parse_out(__doc__))
