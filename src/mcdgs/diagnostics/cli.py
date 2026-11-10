"""The command line of the diagnostics modules, and the files they write.

Each module with a ``main`` runs as ``python -m <module> [--out FOLDER]``. Without ``--out``,
its files go where ``main`` puts them by default, such as the folder of its McStas run.
:func:`save_figure` and :func:`save_text` write a file and say where, and
:func:`markdown_table` lays out the tables the documentation includes.
"""

import argparse
from collections.abc import Iterable, Sequence
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.figure import Figure


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
