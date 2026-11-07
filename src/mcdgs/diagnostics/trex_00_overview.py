"""The whole of T-Rex as McStas sees it: its printout, its instrument file and two views of it.

:func:`main` writes ``print(trex)`` and the McStas instrument file of :class:`.TRex` at the
default control parameters, and draws it with McStasScript's ``show_instrument``, which runs
mcdisplay:

- ``trex_print.txt``: ``print(trex)``, the parameters and the position of every component;
- ``TRex.instr``: the instrument file, written with its build files to ``mcstas_input/TRex/``;
- ``overview_top.png``: the instrument seen from above, z along the beam and x across it;
- ``overview_along_beam.png``: seen from the detector to the source, x across the beam and y up.

Both views have the same scale on their two axes. Run the module for them::

    pixi run python -m mcdgs.diagnostics.trex_00_overview
"""

import shutil
from pathlib import Path
from unittest.mock import patch

import matplotlib.pyplot as plt
from mcstasscript.interface.instr import McStas_instr

from mcdgs.diagnostics.cli import output_folder, parse_out
from mcdgs.instruments.trex.model import TRex
from mcdgs.mcstas.builder import build_mcstas, default_input_dir


def _save_view(path: Path, dpi: int) -> Path:
    """Give the view that show_instrument drew equal axes, save it and close it."""
    plt.gca().set_aspect("equal")
    plt.gcf().savefig(path, dpi=dpi)
    plt.close("all")
    print(f"Saved {path}")
    return path


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
    path = folder / "trex_print.txt"
    path.write_text(f"{trex}\n")
    print(f"Written {path}")
    return path


def draw_views(mcstas_instr: McStas_instr, folder: Path) -> list[Path]:
    """Draw the instrument from above and along the beam, and save both figures.

    show_instrument returns None and has no aspect option, so the aspect is set on the axes it
    drew, and the figure it would show is saved instead.

    Parameters
    ----------
    mcstas_instr:
        The McStasScript instrument.
    folder:
        Folder of the figures.

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
        The paths of the printout, of the instrument file and of the two views.
    """
    trex = TRex()  # the default control parameters
    input_dir = default_input_dir(trex)
    mcstas_instr = build_mcstas(trex, input_path=input_dir)
    mcstas_instr.write_full_instrument()
    instr_file = input_dir / f"{mcstas_instr.name}.instr"

    folder = output_folder(out, input_dir)
    if folder.resolve() != input_dir.resolve():
        instr_file = Path(shutil.copyfile(instr_file, folder / instr_file.name))
    print(f"Written {instr_file}")
    return [write_printout(trex, folder), instr_file, *draw_views(mcstas_instr, folder)]


if __name__ == "__main__":
    main(out=parse_out(__doc__))
