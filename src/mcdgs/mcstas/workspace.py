"""The folders that mcdgs writes into.

A :class:`Workspace` is one folder, its ``root``, with a folder for each kind of file::

    mcdgs_workspace/                 root
        builds/<instrument name>/    the .instr and .c files and the compiled instrument
        runs/<run name>/             what McStas wrote in one run
        results/<name>/              figures, tables and other files made from the runs

It only gives the paths: nothing here knows McStas, an instrument or a run, and no folder is
created before something is written into it.
"""

import os
import shutil
from dataclasses import dataclass
from pathlib import Path

ROOT_VARIABLE = "MCDGS_WORKSPACE"  # the environment variable naming the root
ROOT_NAME = "mcdgs_workspace"  # the root in the current folder, without the variable


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
