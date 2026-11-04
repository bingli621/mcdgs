"""Make the files of the T-Rex page that need McStas, at the default control parameters.

`pixi run docs-mcstas` runs this. Run it again when the instrument changes, and commit what it
writes: the docs build has no McStas, so it shows these as they are.

- docs/instruments/trex/files/TRex.instr: the McStas instrument file, also written with its
  build files to mcstas_input/TRex/ of the repository;
- docs/instruments/trex/figures/overview_top.png: the instrument seen from above, z along the
  beam and x across it;
- docs/instruments/trex/figures/overview_along_beam.png: seen from the detector to the source,
  x across the beam and y up.

The views are McStasScript's show_instrument, which runs mcdisplay.
"""

import shutil
from pathlib import Path
from unittest.mock import patch

import matplotlib

matplotlib.use("Agg")  # no window: the figures only go to files

import matplotlib.pyplot as plt
from mcstasscript.interface.instr import McStas_instr

from mcdgs.instruments.trex.model import TRex
from mcdgs.mcstas.builder import build_mcstas, default_input_dir

REPO = Path(__file__).resolve().parents[1]
DOCS = REPO / "docs" / "instruments" / "trex"
DOCS_INSTR = DOCS / "files" / "TRex.instr"
TOP_VIEW = DOCS / "figures" / "overview_top.png"
BEAM_VIEW = DOCS / "figures" / "overview_along_beam.png"


def write_instr(mcstas_instr: McStas_instr) -> None:
    """Write the .instr file to its input folder, and copy it to the documentation.

    Parameters
    ----------
    mcstas_instr:
        The McStasScript instrument.
    """
    mcstas_instr.write_full_instrument()
    instr_file = Path(mcstas_instr.input_path) / f"{mcstas_instr.name}.instr"
    DOCS_INSTR.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(instr_file, DOCS_INSTR)
    print(f"Written to {instr_file}, copied to {DOCS_INSTR}")


def draw_views(mcstas_instr: McStas_instr) -> None:
    """Draw the instrument from above and along the beam, and save both figures.

    Both have the same scale on their two axes. show_instrument returns None and has no aspect
    option, so the aspect is set on the axes it drew.

    Parameters
    ----------
    mcstas_instr:
        The McStasScript instrument.
    """
    TOP_VIEW.parent.mkdir(parents=True, exist_ok=True)
    # z along the beam, x across it; width and height are in pixels at 100 dpi: 16 x 4 in
    with patch("matplotlib.pyplot.show"):
        mcstas_instr.show_instrument(
            backend="matplotlib_2d", projection="zx", width=1600, height=400
        )
    plt.gca().set_aspect("equal")
    plt.gcf().savefig(TOP_VIEW, dpi=400)  # 6400 px wide, to zoom in on in the documentation
    plt.close("all")
    print(f"Saved {TOP_VIEW}")

    # from the detector to the source; limits in m: [min, max] across, then up
    with patch("matplotlib.pyplot.show"):
        mcstas_instr.show_instrument(
            backend="matplotlib_2d",
            projection="xy",
            limits=[-0.5, 1.5, -1, 1],
            width=600,
            height=600,
        )
    plt.gca().set_aspect("equal")
    plt.gcf().savefig(BEAM_VIEW, dpi=150)
    plt.close("all")
    print(f"Saved {BEAM_VIEW}")


def main() -> None:
    """Write the .instr file and draw both views."""
    trex = TRex()  # the default control parameters
    mcstas_instr = build_mcstas(trex, input_path=default_input_dir(trex, REPO))
    write_instr(mcstas_instr)
    draw_views(mcstas_instr)


if __name__ == "__main__":
    main()
