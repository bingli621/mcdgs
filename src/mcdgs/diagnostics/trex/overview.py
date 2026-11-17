"""The whole of T-Rex as McStas sees it: its printout, its instrument file and two views of it.

:func:`main` writes ``print(trex)``, the McStas instrument file of :class:`.TRex` at the default
control parameters and McStasScript's diagram of how its components are placed, and draws it
with McStasScript's ``show_instrument``, which runs mcdisplay. Each is a function of
:mod:`.diagnostics.drawing`, for any instrument:

- ``trex_print.txt``: ``print(trex)``, the parameters and the position of every component;
- ``TRex.instr``: the instrument file, copied from the build folder of T-Rex;
- ``trex_diagram_legend.png`` and ``trex_diagram.png``: McStasScript's ``show_diagram()``, its
  legend and which component each one is placed and rotated relative to;
- ``overview_top.png``: the instrument seen from above, z along the beam and x across it;
- ``overview_along_beam.png``: seen from the detector to the source, x across the beam and y up,
  with the components from after the source up to before the banana detector, which would
  hide the rest.

Both views have the same scale on their two axes. Run the module for them::

    pixi run python -m mcdgs.diagnostics.trex.overview
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
from mcdgs.instruments.trex.model import Params, TRex
from mcdgs.mcstas import build_mcstas


# --8<-- [start:instrument]
def define_instrument() -> TRex:
    """Return the instrument of this diagnostic: :class:`.TRex` at its default parameters."""
    params = Params()  # the default control parameters
    return TRex(params=params)


# --8<-- [end:instrument]


def main(out: Path | None = None) -> list[Path]:
    """Write the printout and the instrument file of T-Rex, and draw its diagram and both views.

    Parameters
    ----------
    out:
        Folder of the files. If ``None``, ``trex_overview`` in the results folder of
        :meth:`.Workspace.find`.

    Returns
    -------
    :
        The paths of the printout, of the instrument file, of the legend and the diagram, and of
        the two views.
    """
    trex = define_instrument()
    mcstas_instr = build_mcstas(trex)  # in the build folder of T-Rex in the workspace
    folder = output_folder(out, "trex_overview")
    # seen along the beam, the source and the banana detector are wide disks and arcs that
    # cover the guides and the choppers
    first, last = component_range(trex, after="source", before="banana_detector")
    return [
        write_printout(trex, folder / "trex_print.txt"),
        write_instr_file(mcstas_instr, folder),
        *draw_diagram(
            mcstas_instr, folder / "trex_diagram_legend.png", folder / "trex_diagram.png"
        ),
        # from above: z along the beam, x across it, 16 x 4 in; 6400 px wide, to zoom in on in
        # the documentation
        draw_instrument(
            mcstas_instr,
            folder / "overview_top.png",
            dpi=400,
            projection="zx",
            width=1600,
            height=400,
        ),
        # from the detector to the source; limits in m: [min, max] across, then up
        draw_instrument(
            mcstas_instr,
            folder / "overview_along_beam.png",
            dpi=150,
            projection="xy",
            index_min=first,
            index_max=last,
            limits=[-0.5, 1.5, -1, 1],
            width=600,
            height=600,
        ),
    ]


if __name__ == "__main__":
    main(out=parse_out(__doc__))
