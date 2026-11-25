"""The choppers of CSPEC: a picture of each disk, and the acceptance diagram of the cascade.

:func:`main` writes, at the default control parameters:

- ``chopper_<name>.svg``: scippneutron's drawing of each disk, from :func:`.get_choppers`;
- ``cspec_acceptance.png``: the chopper acceptance diagram, from
  :func:`.acceptance_diagram`.

It needs only mcdgs and scippneutron. Run the module for them::

    pixi run python -m mcdgs.diagnostics.cspec.choppers --out docs/instruments/cspec/figures/choppers
"""

from pathlib import Path

from mcdgs.diagnostics import cli
from mcdgs.instruments.cspec.model import CSPEC, Params
from mcdgs.tof import acceptance_diagram, get_choppers


# --8<-- [start:instrument]
def define_instrument() -> CSPEC:
    """Return the instrument of this diagnostic: :class:`.CSPEC` at its default parameters."""
    params = Params()  # the default control parameters
    return CSPEC(params=params)


# --8<-- [end:instrument]


def main(out: Path | None = None) -> list[Path]:
    """Draw every chopper of CSPEC and its acceptance diagram, at the default parameters.

    Parameters
    ----------
    out:
        Folder to write to. If ``None``, ``cspec_choppers`` in the results folder of
        :meth:`.Workspace.find`.

    Returns
    -------
    :
        The paths of the pictures of the disks, in beam order, then of the diagram.
    """
    folder = cli.output_folder(out, "cspec_choppers")
    paths = []

    # --8<-- [start:draw_choppers]
    cspec = define_instrument()
    for name, chopper in get_choppers(cspec).items():
        svg = chopper.make_svg(image_size=300)
        paths.append(cli.save_text(svg, folder / f"chopper_{name}.svg"))
    # --8<-- [end:draw_choppers]

    # --8<-- [start:acceptance]
    figure = acceptance_diagram(cspec, time_max=0.015)
    paths.append(cli.save_figure(figure, folder / "cspec_acceptance.png"))
    # --8<-- [end:acceptance]
    return paths


if __name__ == "__main__":
    main(out=cli.parse_out(__doc__))
