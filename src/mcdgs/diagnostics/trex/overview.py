"""The whole of T-Rex as McStas sees it: its printout, its instrument file and two views of it.

:func:`main` writes ``print(trex)``, the McStas instrument file of :class:`.TRex` at the default
control parameters and McStasScript's diagram of how its components are placed, and draws it
with McStasScript's ``show_instrument``, which runs mcdisplay:

- ``trex_print.txt``: ``print(trex)``, the parameters and the position of every component;
- ``TRex.instr``: the instrument file, written with its build files to ``mcstas_input/TRex/``;
- ``trex_diagram_legend.png`` and ``trex_diagram.png``: McStasScript's ``show_diagram()``, its
  legend and which component each one is placed and rotated relative to;
- ``overview_top.png``: the instrument seen from above, z along the beam and x across it;
- ``overview_along_beam.png``: seen from the detector to the source, x across the beam and y up,
  with the components from after the source up to before the banana detector, which would
  hide the rest.

Both views have the same scale on their two axes. Run the module for them::

    pixi run python -m mcdgs.diagnostics.trex.overview
"""

import shutil
from pathlib import Path
from unittest.mock import patch

import matplotlib.pyplot as plt
from mcstasscript.interface.instr import McStas_instr

from mcdgs.diagnostics.cli import output_folder, parse_out, save_figure, save_text
from mcdgs.instruments.trex.model import TRex
from mcdgs.mcstas.builder import build_mcstas, default_input_dir


def _save_view(path: Path, dpi: int) -> Path:
    """Give the view that show_instrument drew equal axes, save it and close it."""
    plt.gca().set_aspect("equal")
    return save_figure(plt.gcf(), path, dpi=dpi)


def write_printout(trex: TRex, folder: Path) -> Path:
    """Write ``print(trex)`` to ``trex_print.txt``.

    Parameters
    ----------
    trex:
        The instrument.
    folder:
        Folder of the file.

    Returns
    -------
    :
        The path of the file.
    """
    return save_text(f"{trex}\n", folder / "trex_print.txt")


def along_beam_range(trex: TRex) -> tuple[int, int]:
    """Return the indices of the components the view along the beam draws.

    The source and the banana detector are left out: seen along the beam they are wide disks
    and arcs that cover the guides and the choppers.

    Parameters
    ----------
    trex:
        The instrument.

    Returns
    -------
    :
        The index of the first component after the source, and that of the banana detector,
        the first one not drawn.
    """
    names = list(trex.components)
    return names.index("source") + 1, names.index("banana_detector")


def draw_diagram(mcstas_instr: McStas_instr, folder: Path) -> list[Path]:
    """Draw McStasScript's diagram of how the components are placed and its legend, and save them.

    ``show_diagram()`` lists the components in beam order, each with an arrow to the component
    it is placed ``AT ... RELATIVE`` to, in blue, and ``ROTATED ... RELATIVE`` to, in red. It
    draws its legend first, in a figure of its own: the color of each kind of McStas component,
    and the dashed box of a component with a ``WHEN``. It returns None, so the two figures it
    would show are saved instead, without their empty margins.

    Parameters
    ----------
    mcstas_instr:
        The McStasScript instrument.
    folder:
        Folder of the figure.

    Returns
    -------
    :
        The paths of the legend and of the diagram.
    """
    with patch("matplotlib.pyplot.show"):
        mcstas_instr.show_diagram()
    legend, diagram = (plt.figure(number) for number in plt.get_fignums())  # in drawing order
    paths = [folder / "trex_diagram_legend.png", folder / "trex_diagram.png"]
    for figure, path in zip((legend, diagram), paths, strict=True):
        figure.savefig(path, dpi=100, bbox_inches="tight")
        print(f"Saved {path}")
    plt.close("all")
    return paths


def draw_views(
    mcstas_instr: McStas_instr, folder: Path, along_beam_indices: tuple[int, int]
) -> list[Path]:
    """Draw the instrument from above and along the beam, and save both figures.

    show_instrument returns None and has no aspect option, so the aspect is set on the axes it
    drew, and the figure it would show is saved instead.

    Parameters
    ----------
    mcstas_instr:
        The McStasScript instrument.
    folder:
        Folder of the figures.
    along_beam_indices:
        Index of the first component the view along the beam draws, and of the first it does
        not, as :func:`along_beam_range` gives them. The view from above draws them all.

    Returns
    -------
    :
        The paths of the view from above and of the view along the beam.
    """
    # z along the beam, x across it; width and height are in pixels at 100 dpi: 16 x 4 in
    with patch("matplotlib.pyplot.show"):
        mcstas_instr.show_instrument(
            backend="matplotlib_2d", projection="zx", width=1600, height=400
        )
    # 6400 px wide, to zoom in on in the documentation
    top = _save_view(folder / "overview_top.png", dpi=400)

    # from the detector to the source; limits in m: [min, max] across, then up
    with patch("matplotlib.pyplot.show"):
        mcstas_instr.show_instrument(
            backend="matplotlib_2d",
            projection="xy",
            index_min=along_beam_indices[0],
            index_max=along_beam_indices[1],
            limits=[-0.5, 1.5, -1, 1],
            width=600,
            height=600,
        )
    along_beam = _save_view(folder / "overview_along_beam.png", dpi=150)
    return [top, along_beam]


def main(out: Path | None = None) -> list[Path]:
    """Write the printout and the instrument file of T-Rex, and draw both views.

    Parameters
    ----------
    out:
        Folder of the files. If ``None``, ``mcstas_input/TRex/``, where the instrument file is
        written to build it.

    Returns
    -------
    :
        The paths of the printout, of the instrument file, of the legend and the diagram, and of
        the two views.
    """
    trex = TRex()  # the default control parameters
    input_dir = default_input_dir(trex)
    mcstas_instr = build_mcstas(trex, input_path=input_dir)
    mcstas_instr.write_full_instrument()
    instr_file = input_dir / f"{mcstas_instr.name}.instr"

    folder = output_folder(out, input_dir)
    if folder.resolve() != input_dir.resolve():
        instr_file = Path(shutil.copyfile(instr_file, folder / instr_file.name))
    print(f"Saved {instr_file}")
    return [
        write_printout(trex, folder),
        instr_file,
        *draw_diagram(mcstas_instr, folder),
        *draw_views(mcstas_instr, folder, along_beam_range(trex)),
    ]


if __name__ == "__main__":
    main(out=parse_out(__doc__))
