"""The command line of the diagnostics modules: where their files go.

Each module with a ``main`` runs as ``python -m <module> [--out FOLDER]``. Without ``--out``,
its files go where ``main`` puts them by default, such as the folder of its McStas run.
"""

import argparse
from collections.abc import Sequence
from pathlib import Path


def parse_out(description: str | None, argv: Sequence[str] | None = None) -> Path | None:
    """Read ``--out`` from the command line.

    Parameters
    ----------
    description:
        What the module does, for ``--help``: its docstring.
    argv:
        The arguments. If ``None``, those of the command line.

    Returns
    -------
    :
        The folder given with ``--out``, or ``None``.
    """
    first_line = description.strip().splitlines()[0] if description else None
    parser = argparse.ArgumentParser(description=first_line)
    parser.add_argument("--out", type=Path, help="folder to write the files to")
    return parser.parse_args(argv).out


def output_folder(out: Path | None, default: Path) -> Path:
    """Return the folder to write to, created with its parents.

    Parameters
    ----------
    out:
        The folder given with ``--out``, or ``None``.
    default:
        The folder to use without ``out``.

    Returns
    -------
    :
        ``out``, or ``default`` without it.
    """
    folder = default if out is None else out
    folder.mkdir(parents=True, exist_ok=True)
    return folder
