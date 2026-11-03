"""Write the .instr file of T-Rex, with its default control parameters, to mcstas_input/TRex/ of
the repository, print it, and plot the instrument with McStas's mcdisplay: seen from above, and
along the beam. The .instr file and both figures are saved for the documentation.
"""

import shutil
from pathlib import Path

import matplotlib.pyplot as plt

from mcdgs.instruments.trex.model import TRex
from mcdgs.mcstas.builder import build_mcstas, default_input_dir

WORK_DIR = Path(__file__).resolve().parents[2]  # the repository, wherever the script runs from
DOCS = WORK_DIR / "docs" / "instruments" / "trex"
DOCS_INSTR = DOCS / "files" / "TRex.instr"  # the copy the documentation shows
FIGURES = DOCS / "figures"
TOP_VIEW = FIGURES / "overview_top.png"
BEAM_VIEW = FIGURES / "overview_along_beam.png"


if __name__ == "__main__":
    trex = TRex()
    mcstas_instr = build_mcstas(trex, input_path=default_input_dir(trex, WORK_DIR))
    mcstas_instr.write_full_instrument()
    instr_file = Path(mcstas_instr.input_path) / f"{trex.name}.instr"
    print(instr_file.read_text())
    print(f"Written to {instr_file}")
    DOCS_INSTR.parent.mkdir(exist_ok=True)
    shutil.copyfile(instr_file, DOCS_INSTR)
    print(f"Copied to {DOCS_INSTR}")

    # show_instrument calls plt.show() itself, which in a script blocks on each figure: in
    # interactive mode it returns at once, and the plt.show() at the end shows them all
    plt.ion()
    FIGURES.mkdir(exist_ok=True)
    # z along the beam, x across it; limits are [zmin, zmax, xmin, xmax] in m
    # width and height are in pixels at 100 dpi: a figure of 16 x 4 in
    mcstas_instr.show_instrument(backend="matplotlib_2d", projection="zx", width=1600, height=400)
    # show_instrument returns None and has no aspect option: set it on the axes it drew
    plt.gca().set_aspect("equal")
    plt.gcf().savefig(TOP_VIEW, dpi=400)  # 6400 px wide, to zoom in on in the documentation
    print(f"Saved {TOP_VIEW}")
    mcstas_instr.show_instrument(
        backend="matplotlib_2d", projection="xy", limits=[-0.5, 1.5, -1, 1], width=600, height=600
    )  # along the beam
    plt.gca().set_aspect("equal")
    plt.gcf().savefig(BEAM_VIEW, dpi=150)
    print(f"Saved {BEAM_VIEW}")
    plt.ioff()
    plt.show()
