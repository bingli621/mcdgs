"""Run folders that a :class:`.Simulation` writes into.

Each run writes into its own folder in the output folder (:func:`default_output_dir`). By
default the folder is named after the start time and what the run does, e.g.
``20261022-143005_scan_omega``.

- :func:`run_folder` picks that folder and clears the way for the run to write it.
- :func:`list_runs`, :func:`latest_run` and :func:`latest_beam` find finished runs.
- :func:`clean_runs` removes older runs, never a ``run_to_sample`` one.
"""

import re
import shutil
from datetime import datetime
from pathlib import Path


def default_output_dir(work_dir: str | Path | None = None) -> Path:
    """Return the folder of the simulation output, McStasScript's ``output_path``.

    It holds one folder per run, and each run creates it if missing.

    Parameters
    ----------
    work_dir:
        Work folder. If ``None``, the current folder is used.

    Returns
    -------
    :
        ``work_dir/mcstas_output``.
    """
    return Path(work_dir or Path.cwd()) / "mcstas_output"


def _check_name(folder_name: object) -> None:
    """Check that ``folder_name`` is a plain folder name.

    A plain name has letters, digits, ``.``, ``_`` and ``-``, and starts with a letter or digit.
    A run removes its folder before writing it, so a name like ``..``, ``/home`` or ``a/b`` must
    never get through: it would remove a folder outside the output folder.

    Raises
    ------
    ValueError
        If ``folder_name`` is not a string or not a plain folder name.
    """
    safe = isinstance(folder_name, str) and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", folder_name)
    if not safe:
        raise ValueError(
            f"folder_name {folder_name!r}: one folder name of letters, digits, '.', '_' and '-'"
        )


def _new_name(output_dir: Path, mode: str) -> str:
    """Return a new run folder name, not yet taken in ``output_dir``.

    The name is the local time now and ``mode``, e.g. ``20261022-143005_scan_omega``. A run
    started in the same second as an earlier one gets ``-2``, ``-3``, ... appended, so it never
    overwrites the earlier run.
    """
    name = unique = f"{datetime.now().astimezone():%Y%m%d-%H%M%S}_{mode}"
    n = 1
    while (output_dir / unique).exists():
        n += 1
        unique = f"{name}-{n}"
    return unique


def run_folder(output_dir: Path, folder_name: str | None, mode: str) -> Path:
    """Return the folder a run writes into, and clear the way for it.

    The run folder itself is not created, since McStas refuses to write into an existing folder.
    ``output_dir`` is created with its parents, since McStas does not create them.

    Parameters
    ----------
    output_dir:
        Output folder.
    folder_name:
        Name of the run folder. A folder of that name from an earlier run is removed, so the run
        replaces it. If ``None``, a new name from the time now and ``mode`` is used, so each run
        keeps its own output.
    mode:
        What the run does, e.g. ``"run"``, ``"scan_omega"``, ``"run_to_sample"``.

    Returns
    -------
    :
        ``output_dir / folder_name``.

    Raises
    ------
    ValueError
        If ``folder_name`` is not a plain folder name of letters, digits, ``.``, ``_`` and ``-``.
    """
    if folder_name is None:
        folder_name = _new_name(output_dir, mode)
    _check_name(folder_name)
    path = output_dir / folder_name
    shutil.rmtree(path, ignore_errors=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    return path


def list_runs(output_dir: str | Path | None = None) -> list[str]:
    """List the finished run folders, newest first.

    A run is finished when its folder has a ``mccode.h5``.

    Parameters
    ----------
    output_dir:
        Output folder. If ``None``, :func:`default_output_dir` is used.

    Returns
    -------
    :
        Folder names, e.g. ``["20261022-143005_scan_omega", "20261022-120000_run_to_sample"]``.
        Empty if ``output_dir`` does not exist.
    """
    output_dir = Path(output_dir or default_output_dir())
    if not output_dir.is_dir():
        return []
    runs = [d for d in output_dir.iterdir() if (d / "mccode.h5").is_file()]
    runs.sort(key=lambda d: (d / "mccode.h5").stat().st_mtime, reverse=True)
    return [run.name for run in runs]


def latest_run(output_dir: str | Path | None = None, pattern: str = "*") -> Path:
    """Return the newest finished run folder whose name matches a pattern.

    Parameters
    ----------
    output_dir:
        Output folder. If ``None``, :func:`default_output_dir` is used.
    pattern:
        Glob pattern for the folder name, e.g. ``"*_scan_omega"``, ``"*_run_to_sample"``,
        ``"*_from_sample"``.

    Returns
    -------
    :
        The run folder.

    Raises
    ------
    FileNotFoundError
        If no finished run matches.
    """
    output_dir = Path(output_dir or default_output_dir())
    for name in list_runs(output_dir):
        if Path(name).match(pattern):
            return output_dir / name
    raise FileNotFoundError(f"no finished run {pattern!r} in {output_dir}")


def latest_beam(output_dir: str | Path | None = None) -> Path | None:
    """Return the newest beam file of the finished runs.

    Parameters
    ----------
    output_dir:
        Output folder. If ``None``, :func:`default_output_dir` is used.

    Returns
    -------
    :
        The ``beam.mcpl.gz`` of the newest run that has one (a ``run_to_sample``), or ``None``.
    """
    output_dir = Path(output_dir or default_output_dir())
    beams = [output_dir / run / "beam.mcpl.gz" for run in list_runs(output_dir)]
    return next((beam for beam in beams if beam.is_file()), None)


def monitor_file(folder: str | Path, output_dir: str | Path | None = None) -> Path | None:
    """Return the file with the incident monitor that normalizes a run.

    A run from the sample (a folder name ending in ``_from_sample``) has no monitor in its own
    file, so it is normalized by the newest ``run_to_sample``.

    Parameters
    ----------
    folder:
        The run folder.
    output_dir:
        Where to look for the ``run_to_sample``. If ``None``, ``folder``'s parent.

    Returns
    -------
    :
        For a run from the sample, the ``mccode.h5`` of the newest ``run_to_sample``. Otherwise
        ``None``, as the run has its own monitor.

    Raises
    ------
    FileNotFoundError
        If a run from the sample finds no ``run_to_sample``.
    """
    folder = Path(folder)
    if not folder.name.endswith("_from_sample"):
        return None
    return latest_run(output_dir or folder.parent, "*_run_to_sample") / "mccode.h5"


def clean_runs(output_dir: str | Path | None = None, keep: int = 5) -> list[Path]:
    """Remove all but the newest finished runs.

    Never removed, whatever their age:

    - a ``run_to_sample`` run, i.e. any folder with a ``beam.mcpl.gz``. Its beam file is the
      expensive simulation every ``run_from_sample`` and ``scan_from_sample`` starts from, and
      its ``mccode.h5`` has the incident monitor that normalizes them. Delete one by hand when
      it is no longer needed.
    - a folder without a ``mccode.h5``: an unfinished run, or anything else in ``output_dir``.

    Parameters
    ----------
    output_dir:
        Output folder. If ``None``, :func:`default_output_dir` is used.
    keep:
        Number of newest finished runs to keep. ``run_to_sample`` runs count toward it, so
        fewer other runs may be kept.

    Returns
    -------
    :
        The removed folders.
    """
    output_dir = Path(output_dir or default_output_dir())
    runs = [output_dir / run for run in list_runs(output_dir)]
    old = [run for run in runs[keep:] if not (run / "beam.mcpl.gz").is_file()]
    for run in old:
        shutil.rmtree(run)
    return old
