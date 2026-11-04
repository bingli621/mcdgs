"""What the scripts that make the files of the documentation share.

- docs_build.py makes what needs only mcdgs and tof, for every build of the site;
- docs_mcstas_*.py make what needs McStas: run by hand, and their files committed.

The pixi tasks run them with MPLBACKEND=Agg: the figures only go to files.
"""

from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
TREX_DOCS = REPO / "docs" / "instruments" / "trex"
SEED = 1  # the same neutrons in every run, so that the pages and their files agree
