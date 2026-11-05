"""Make the figure of the T-Rex source page, which needs McStas.

`pixi run docs-mcstas` runs this with docs_mcstas_overview.py and docs_mcstas_bender.py. Run
it when the source changes, and commit what it writes: the docs build has no McStas, so it shows
the figure as it is.

- docs/instruments/trex/figures/trex_source.png: McStas and tof at the entrance of guide 1, from
  ``mcdgs.diagnostics.trex.source.main`` with its default numbers of neutrons.
"""

import shutil

from docs_common import SEED, TREX_DOCS

from mcdgs.diagnostics.trex.source import main as compare_source

FIGURE = TREX_DOCS / "figures" / "trex_source.png"


def main() -> None:
    """Run both models of the source and copy their figure to the documentation."""
    figure = compare_source(seed=SEED)
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(figure, FIGURE)
    print(f"Copied to {FIGURE}")


if __name__ == "__main__":
    main()
