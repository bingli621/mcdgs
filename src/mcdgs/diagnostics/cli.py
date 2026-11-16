"""The command line of the diagnostics modules, and the files they write.

Each module with a ``main`` runs as ``python -m <module> [--out FOLDER]``. Without ``--out``,
its files go into a folder of its own in the results folder of the workspace
(:func:`output_folder`), apart from the McStas runs they are made from.
:func:`save_figure` and :func:`save_text` write a file and say where, and
:func:`markdown_table` lays out the tables the documentation includes.
:func:`remove_mcdisplay_folders` clears what drawing an instrument leaves in its build folder.
"""

import argparse
import shutil
from collections.abc import Iterable, Sequence
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from mcdgs.mcstas import Workspace


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


def output_folder(out: Path | None, name: str) -> Path:
    """Return the folder to write to, created with its parents.

    Parameters
    ----------
    out:
        The folder given with ``--out``, or ``None``.
    name:
        The name of the folder to use without ``out``, in the results folder of the workspace:
        that of the module, e.g. ``"trex_choppers"``.

    Returns
    -------
    :
        ``out``, or without it ``name`` in the results folder of :meth:`.Workspace.find`.
    """
    folder = Workspace.find().result_folder(name) if out is None else out
    folder.mkdir(parents=True, exist_ok=True)
    return folder


def remove_mcdisplay_folders(build_folder: str | Path) -> list[Path]:
    """Remove the folders mcdisplay left in the build folder of an instrument.

    McStasScript's ``show_instrument`` runs mcdisplay, which writes a new folder
    ``<instrument>_mcdisplay``, ``<instrument>_mcdisplay_0``, ... beside the instrument file at
    every call. They are not needed once the instrument is drawn.

    Parameters
    ----------
    build_folder:
        The build folder of the instrument, McStasScript's ``input_path``.

    Returns
    -------
    :
        The removed folders.
    """
    folders = sorted(p for p in Path(build_folder).glob("*_mcdisplay*") if p.is_dir())
    for folder in folders:
        shutil.rmtree(folder)
    return folders


def save_figure(fig: Figure, path: Path, dpi: int = 150) -> Path:
    """Save a figure, close it and say where it went.

    Parameters
    ----------
    fig:
        The figure.
    path:
        The file, whose extension gives the format.
    dpi:
        Its resolution, in dots per inch.

    Returns
    -------
    :
        ``path``.
    """
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    print(f"Saved {path}")
    return path


def save_text(text: str, path: Path) -> Path:
    """Write a text file and say where it went.

    Parameters
    ----------
    text:
        The text.
    path:
        The file, overwritten if it exists.

    Returns
    -------
    :
        ``path``.
    """
    path.write_text(text)
    print(f"Saved {path}")
    return path


def markdown_table(header: Sequence[str], rows: Iterable[Sequence[object]]) -> str:
    """Return a Markdown table, each value as it is formatted.

    A ``-`` in a value becomes a minus sign, ``−``, as in the documentation, and an empty value
    an empty cell.

    Parameters
    ----------
    header:
        The name of each column.
    rows:
        The values of each row, one per column, already formatted, such as ``f"{x:.3f} ms"``.

    Returns
    -------
    :
        The table, one line per row, ending in a newline.

    Examples
    --------
    >>> print(markdown_table(["Chopper", "Frequency"], [["`ps2`", "-189 Hz"]]), end="")
    | Chopper | Frequency |
    | ------- | --------- |
    | `ps2` | −189 Hz |
    """
    lines = [header, ["-" * len(name) for name in header]]
    lines += [[str(value).replace("-", "\u2212") for value in row] for row in rows]
    return "".join("|" + "|".join(f" {v} " if v else " " for v in line) + "|\n" for line in lines)
