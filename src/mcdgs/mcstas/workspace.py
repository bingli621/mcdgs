"""The folders that mcdgs writes into.

A :class:`Workspace` is one folder, its ``root``, with a folder for each kind of file::

    mcdgs_workspace/                 root
        builds/<instrument name>/    the .instr and .c files and the compiled instrument
        runs/<run name>/             what McStas wrote in one run
        results/<name>/              figures, tables and other files made from the runs

It gives the paths, and looks after what is in them:

- :meth:`Workspace.new_run_folder` picks the folder of a run and clears the way for the run to
  write it. By default the folder is named after the start time and what the run does, e.g.
  ``20261022-143005_scan_omega``.
- :meth:`Workspace.list_runs` and :meth:`Workspace.latest_run` find the finished runs, by their
  records (:class:`.Run`).
- :meth:`Workspace.clean_runs` removes older runs, and :meth:`Workspace.clean_builds` the builds.

Nothing here knows McStas or an instrument, and no folder is created before something is written
into it.
"""

import os
import re
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from mcdgs.mcstas.runs import RUN_FILE, Run, RunType

ROOT_VARIABLE = "MCDGS_WORKSPACE"  # the environment variable naming the root
ROOT_NAME = "mcdgs_workspace"  # the root in the current folder, without the variable


def _check_name(folder_name: object) -> None:
    """Check that ``folder_name`` is a plain folder name.

    A plain name has letters, digits, ``.``, ``_`` and ``-``, and starts with a letter or digit.
    A run removes its folder before writing it, so a name like ``..``, ``/home`` or ``a/b`` must
    never get through: it would remove a folder outside the runs folder.

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


def _new_name(runs_folder: Path, label: str) -> str:
    """Return a new run folder name, not yet taken in ``runs_folder``.

    The name is the local time now and ``label``, e.g. ``20261022-143005_scan_omega``. A run
    started in the same second as an earlier one gets ``-2``, ``-3``, ... appended, so it never
    overwrites the earlier run.
    """
    name = unique = f"{datetime.now().astimezone():%Y%m%d-%H%M%S}_{label}"
    n = 1
    while (runs_folder / unique).exists():
        n += 1
        unique = f"{name}-{n}"
    return unique


@dataclass(frozen=True)
class Workspace:
    """The folders of the builds, the runs and the results, in one root folder.

    Parameters
    ----------
    root:
        The folder that holds them.
    """

    root: Path

    def __post_init__(self) -> None:
        object.__setattr__(self, "root", Path(self.root))

    @classmethod
    def find(cls, root: str | Path | None = None) -> "Workspace":
        """Return the workspace to use.

        Parameters
        ----------
        root:
            Its root folder. If ``None``, the folder the environment variable
            ``MCDGS_WORKSPACE`` names, else ``mcdgs_workspace`` in the current folder.

        Returns
        -------
        :
            The workspace. Its root is fixed here, so a later change of the current folder
            does not move it.
        """
        if root is None:
            root = os.environ.get(ROOT_VARIABLE) or Path.cwd() / ROOT_NAME
        return cls(Path(root).expanduser())

    @property
    def builds_folder(self) -> Path:
        """The folder of the builds, ``root/builds``."""
        return self.root / "builds"

    @property
    def runs_folder(self) -> Path:
        """The folder of the runs, ``root/runs``."""
        return self.root / "runs"

    @property
    def results_folder(self) -> Path:
        """The folder of the results, ``root/results``."""
        return self.root / "results"

    def build_folder(self, instrument_name: str) -> Path:
        """Return the folder of the McStas-generated files of an instrument.

        This is McStasScript's ``input_path``, where mcrun runs, and holds the .instr, .c and
        compiled binary. Each instrument name has its own folder, so two instruments do not
        overwrite each other.

        Parameters
        ----------
        instrument_name:
            The name of the instrument.

        Returns
        -------
        :
            ``builds_folder/<instrument name>``.
        """
        return self.builds_folder / instrument_name

    def run_folder(self, run_name: str) -> Path:
        """Return the folder of a run.

        Parameters
        ----------
        run_name:
            The name of the run folder.

        Returns
        -------
        :
            ``runs_folder/<run name>``.
        """
        return self.runs_folder / run_name

    def result_folder(self, name: str) -> Path:
        """Return a folder for files made from the runs, such as figures and tables.

        Parameters
        ----------
        name:
            The name of the folder, e.g. that of the diagnostic that writes it.

        Returns
        -------
        :
            ``results_folder/<name>``.
        """
        return self.results_folder / name

    def new_run_folder(self, folder_name: str | None, label: str) -> Path:
        """Return the folder a run writes into, and clear the way for it.

        The run folder itself is not created, since McStas refuses to write into an existing
        folder. ``runs_folder`` is created with its parents, since McStas does not create them.

        Parameters
        ----------
        folder_name:
            Name of the run folder. A folder of that name from an earlier run is removed, so the
            run replaces it. If ``None``, a new name from the time now and ``label`` is used, so
            each run keeps its own output.
        label:
            What the run does, for the new name, e.g. ``"run"``, ``"scan_omega"``,
            ``"run_to_sample"``.

        Returns
        -------
        :
            ``runs_folder/<folder name>``.

        Raises
        ------
        ValueError
            If ``folder_name`` is not a plain folder name of letters, digits, ``.``, ``_`` and
            ``-``.
        """
        if folder_name is None:
            folder_name = _new_name(self.runs_folder, label)
        _check_name(folder_name)
        path = self.run_folder(folder_name)
        shutil.rmtree(path, ignore_errors=True)
        self.runs_folder.mkdir(parents=True, exist_ok=True)
        return path

    def list_runs(self) -> list[Run]:
        """List the finished runs, newest first.

        A run is finished when its folder has a ``run.json``. Any other folder is left out: an
        unfinished run, or something else in ``runs_folder``.

        Returns
        -------
        :
            The runs, by the time they started, the newest first. Empty if there is no
            ``runs_folder``.
        """
        if not self.runs_folder.is_dir():
            return []
        runs = [Run.load(d) for d in self.runs_folder.iterdir() if (d / RUN_FILE).is_file()]
        return sorted(runs, key=lambda run: run.started, reverse=True)

    def latest_run(self, run_type: RunType | None = None, instrument: str | None = None) -> Run:
        """Return the newest finished run of a type and of an instrument.

        Parameters
        ----------
        run_type:
            What the run did. If ``None``, any type.
        instrument:
            The name of the instrument. If ``None``, any instrument.

        Returns
        -------
        :
            The run.

        Raises
        ------
        FileNotFoundError
            If no finished run matches.
        """
        for run in self.list_runs():
            if run_type in (None, run.run_type) and instrument in (None, run.instrument):
                return run
        wanted = [
            f"{what} {value}"
            for what, value in (("type", run_type), ("instrument", instrument))
            if value is not None
        ]
        of = f" of {' and '.join(wanted)}" if wanted else ""
        raise FileNotFoundError(f"no finished run{of} in {self.runs_folder}")

    def clean_runs(self, keep: int = 5, run_type: RunType | None = None) -> list[Path]:
        """Remove all but the newest finished runs.

        Only a finished run is ever removed: a folder without a ``run.json``, an unfinished run
        or anything else in ``runs_folder``, is left as it is.

        Parameters
        ----------
        keep:
            Number of newest runs to keep.
        run_type:
            The type of the runs to clean: the newest ``keep`` of that type are kept, its older
            ones removed, and the runs of every other type left as they are. If ``None``, the
            runs of all types are cleaned together, but a ``run_to_sample`` is never removed:
            its beam file is the expensive simulation every run from the sample starts from, and
            its ``mccode.h5`` has the incident monitor that normalizes them. It still counts
            toward ``keep``, so fewer other runs may be kept. Clean those with
            ``run_type=RunType.RUN_TO_SAMPLE``.

        Returns
        -------
        :
            The removed folders, the newest first.
        """
        runs = self.list_runs()
        if run_type is None:
            old = [run for run in runs[keep:] if run.run_type is not RunType.RUN_TO_SAMPLE]
        else:
            old = [run for run in runs if run.run_type is run_type][keep:]
        for run in old:
            shutil.rmtree(run.folder)
        return [run.folder for run in old]

    def clean_builds(self) -> list[Path]:
        """Remove the build folder of every instrument.

        Nothing is lost: every run builds its instrument anew.

        Returns
        -------
        :
            The removed folders. Empty if there is no ``builds_folder``.
        """
        if not self.builds_folder.is_dir():
            return []
        builds = sorted(path for path in self.builds_folder.iterdir() if path.is_dir())
        for build in builds:
            shutil.rmtree(build)
        return builds
