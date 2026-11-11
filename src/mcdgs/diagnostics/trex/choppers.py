"""The choppers of T-Rex: a picture of each disk, and the acceptance diagram of the cascade.

:func:`main` writes, at the default control parameters:

- ``chopper_<name>.svg``: scippneutron's drawing of each disk, from :func:`.get_choppers`;
- ``trex_acceptance.png``: the chopper acceptance diagram, from
  :func:`.acceptance_diagram`.

It needs only mcdgs and scippneutron. Run the module for them::

    pixi run python -m mcdgs.diagnostics.trex.choppers --out docs/instruments/trex/figures/choppers
"""

from pathlib import Path

from mcdgs.diagnostics.cli import output_folder, parse_out, save_figure, save_text
from mcdgs.instruments.trex.model import TRex
from mcdgs.mcstas.run_folders import default_output_dir
from mcdgs.tof.acceptance import acceptance_diagram
from mcdgs.tof.builder import get_choppers


def main(out: Path | None = None) -> list[Path]:
    """Draw every chopper of T-Rex and its acceptance diagram, at the default parameters.

    Parameters
    ----------
    out:
        Folder to write to. If ``None``, ``mcstas_output/trex_choppers``, beside the McStas runs.

    Returns
    -------
    :
        The paths of the pictures of the disks, in beam order, then of the diagram.
    """
    folder = output_folder(out, default_output_dir() / "trex_choppers")
    paths = []

    # --8<-- [start:draw_choppers]
    trex = TRex()  # the default control parameters
    for name, chopper in get_choppers(trex).items():
        svg = chopper.make_svg(image_size=300)
        paths.append(save_text(svg, folder / f"chopper_{name}.svg"))
    # --8<-- [end:draw_choppers]

    # --8<-- [start:acceptance]
    figure = acceptance_diagram(trex, time_max=0.015)
    paths.append(save_figure(figure, folder / "trex_acceptance.png"))
    # --8<-- [end:acceptance]
    return paths


if __name__ == "__main__":
    main(out=parse_out(__doc__))
