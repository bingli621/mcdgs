"""The choppers of T-Rex: a picture of each disk, and the acceptance diagram of the cascade.

:func:`main` writes, at the default control parameters:

- ``chopper_<name>.svg``: scippneutron's drawing of each disk, from :func:`.get_choppers`;
- ``trex_acceptance.png``: the chopper acceptance diagram, from
  :meth:`.TRex.acceptance_diagram`.

It needs only mcdgs and scippneutron. Run the module for them::

    pixi run python -m mcdgs.diagnostics.trex_05_choppers --out docs/instruments/trex/figures
"""

from pathlib import Path

import matplotlib.pyplot as plt

from mcdgs.diagnostics.cli import output_folder, parse_out
from mcdgs.instruments.trex.model import TRex
from mcdgs.tof.builder import get_choppers


def main(out: Path | None = None) -> list[Path]:
    """Draw every chopper of T-Rex and its acceptance diagram, at the default parameters.

    Parameters
    ----------
    out:
        Folder to write to. If ``None``, the current folder.

    Returns
    -------
    :
        The paths of the pictures of the disks, in beam order, then of the diagram.
    """
    folder = output_folder(out, Path.cwd())
    paths = []

    # --8<-- [start:draw_choppers]
    trex = TRex()  # the default control parameters
    for name, chopper in get_choppers(trex).items():
        path = folder / f"chopper_{name}.svg"
        path.write_text(chopper.make_svg(image_size=300))
        # --8<-- [end:draw_choppers]
        paths.append(path)

    # --8<-- [start:acceptance]
    figure = trex.acceptance_diagram()
    figure.savefig(folder / "trex_acceptance.png", dpi=150)
    # --8<-- [end:acceptance]
    plt.close(figure)
    paths.append(folder / "trex_acceptance.png")

    print(f"Written to {folder}")
    return paths


if __name__ == "__main__":
    main(out=parse_out(__doc__))
