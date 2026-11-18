"""The whole of CSPEC as McStas sees it: its printout, its instrument file and two views of it.

:func:`main` writes ``print(cspec)``, the McStas instrument file of :class:`.CSPEC` at the
default control parameters and McStasScript's diagram of how its components are placed, and
draws it with McStasScript's ``show_instrument``, which runs mcdisplay. Each is a function of
:mod:`.diagnostics.drawing`, for any instrument:

- ``cspec_print.txt``: ``print(cspec)``, the parameters and the position of every component;
- ``CSPEC.instr``: the instrument file, copied from the build folder of CSPEC;
- ``cspec_diagram_legend.png`` and ``cspec_diagram.png``: McStasScript's ``show_diagram()``, its
  legend and which component each one is placed and rotated relative to;
- ``overview_top.png``: the instrument seen from above, z along the beam and x across it,
  with x stretched;
- ``overview_side.png``: seen from the side, z along the beam and y up, with the drop of the
  S-bender;
- ``overview_along_beam.png``: seen from the sample to the source, x across the beam and y up,
  with the components from after the source to ``sample_pos``.

The view along the beam has the same scale on its two axes. Run the module for them::

    pixi run python -m mcdgs.diagnostics.cspec.overview
"""

from pathlib import Path

from mcdgs.diagnostics.cli import output_folder, parse_out
from mcdgs.diagnostics.drawing import (
    component_range,
    draw_diagram,
    draw_instrument,
    write_instr_file,
    write_printout,
)
from mcdgs.instruments.cspec.model import CSPEC, Params
from mcdgs.mcstas import build_mcstas


# --8<-- [start:instrument]
def define_instrument() -> CSPEC:
    """Return the instrument of this diagnostic: :class:`.CSPEC` at its default parameters."""
    params = Params()  # the default control parameters
    return CSPEC(params=params)


# --8<-- [end:instrument]


def main(out: Path | None = None) -> list[Path]:
    """Write the printout and the instrument file of CSPEC, and draw its diagram and its views.

    Parameters
    ----------
    out:
        Folder of the files. If ``None``, ``cspec_overview`` in the results folder of
        :meth:`.Workspace.find`.

    Returns
    -------
    :
        The paths of the printout, of the instrument file, of the legend and the diagram, and of
        the three views.
    """
    cspec = define_instrument()
    mcstas_instr = build_mcstas(cspec)  # in the build folder of CSPEC in the workspace
    folder = output_folder(out, "cspec_overview")
    # seen along the beam, the source is a wide disk that covers the guides and the choppers
    first, last = component_range(cspec, after="source", before="sample_pos")
    return [
        write_printout(cspec, folder / "cspec_print.txt"),
        write_instr_file(mcstas_instr, folder),
        *draw_diagram(
            mcstas_instr, folder / "cspec_diagram_legend.png", folder / "cspec_diagram.png"
        ),
        # from above: z along the beam, x across it, 16 x 4 in; 6400 px wide, to zoom in on in
        # the documentation. x stretched, or the 0.35 m the curve takes the beam across would
        # not show
        draw_instrument(
            mcstas_instr,
            folder / "overview_top.png",
            dpi=400,
            equal_aspect=False,
            projection="zx",
            limits=[-5.0, 165.0, -0.2, 0.45],  # along, then across
            width=1600,
            height=400,
        ),
        # from the side: z along the beam, y up, likewise; y stretched, or the drop of 0.42 m
        # over 160 m would not show
        draw_instrument(
            mcstas_instr,
            folder / "overview_side.png",
            dpi=400,
            equal_aspect=False,
            projection="zy",
            limits=[-5.0, 165.0, -0.6, 0.2],  # along, then up
            width=1600,
            height=400,
        ),
        # from the sample to the source; limits in m: [min, max] across, then up
        draw_instrument(
            mcstas_instr,
            folder / "overview_along_beam.png",
            dpi=150,
            projection="xy",
            index_min=first,
            index_max=last,
            limits=[-0.6, 1.0, -1.3, 0.8],
            width=600,
            height=600,
        ),
    ]


if __name__ == "__main__":
    main(out=parse_out(__doc__))
