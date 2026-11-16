"""The runs of a :class:`.Simulation`: their folders and their records.

Each run writes into its own folder in the runs folder (:attr:`.Workspace.runs_folder`). By
default the folder is named after the start time and what the run does, e.g.
``20261022-143005_scan_omega``.

A :class:`Run` says what one run did: which instrument, with which :class:`RunType` and which
values of the instrument parameters, and when it started and finished. It is saved as
``run.json`` in the run folder once McStas has finished (:meth:`Run.finish`), so a folder with
a ``run.json`` is a finished run, and it is read back from there, also in a later session and
without the instrument.

- :func:`recorded_instr_params` gives the instrument parameters of a run as its record keeps
  them.
- :func:`new_run_folder` picks the folder of a run and clears the way for the run to write it.
- :func:`list_runs` and :func:`latest_run` find finished runs, and :func:`monitor_file` the
  file a run is normalized by.
- :func:`clean_runs` removes older runs.

Nothing here knows McStas, an :class:`.Instrument` or a :class:`.Workspace`: the runs are found
by their folder, and what a run did is read from its record, never from the name of its folder.
"""

import json
import re
import shutil
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timedelta
from enum import StrEnum
from numbers import Real
from pathlib import Path

from mcdgs.core.params import InstrParam

RUN_FILE = "run.json"  # the record of a run, in its folder
DATA_FILE = "mccode.h5"  # the McStas NeXus file of a run
BEAM_FILE = "beam.mcpl.gz"  # the beam a run to the sample dumps


class RunType(StrEnum):
    """What a run does: the method of :class:`.Simulation` that made it."""

    RUN = "run"
    SCAN = "scan"
    RUN_TO_SAMPLE = "run_to_sample"
    RUN_FROM_SAMPLE = "run_from_sample"
    SCAN_FROM_SAMPLE = "scan_from_sample"


@dataclass(frozen=True)
class Run:
    """One finished run and what it simulated.

    Parameters
    ----------
    folder:
        The run folder. It is where the record was found, and is not saved in it, so a run
        folder can be moved or renamed.
    instrument:
        The name of the instrument.
    run_type:
        What the run did.
    started:
        When the run started, with its time zone.
    instr_params:
        The value of every instrument parameter in the run, by McStas name. Of a scan, the
        ``variable`` is left out: its values are in ``values``.
    finished:
        When McStas finished, with its time zone. ``None`` while the run has not finished, and
        in a record saved before the finish was recorded.
    ncount:
        The number of neutrons, of a scan per value. ``None`` from the sample, where it is the
        number in the beam file.
    seed:
        The random seed given. ``None`` if McStas picked one.
    variable:
        The instrument parameter of a scan, else ``None``.
    values:
        The values of ``variable``, in the order of the entries of the file, else ``None``.
    beam_run:
        The name of the run to the sample whose beam a run from the sample started from, in
        the same runs folder, else ``None``.
    """

    folder: Path
    instrument: str
    run_type: RunType
    started: datetime
    instr_params: dict[str, float | str]
    finished: datetime | None = None
    ncount: float | None = None
    seed: int | None = None
    variable: str | None = None
    values: list[float] | None = None
    beam_run: str | None = None

    @property
    def name(self) -> str:
        """The name of the run: that of its folder."""
        return self.folder.name

    @property
    def data_file(self) -> Path:
        """The McStas NeXus file of the run, ``folder/mccode.h5``."""
        return self.folder / DATA_FILE

    @property
    def beam_file(self) -> Path | None:
        """The beam file of a run to the sample, ``folder/beam.mcpl.gz``, else ``None``."""
        return self.folder / BEAM_FILE if self.run_type is RunType.RUN_TO_SAMPLE else None

    @property
    def duration(self) -> timedelta | None:
        """The time the run took, from ``started`` to ``finished``: ``None`` without ``finished``.

        It is the time of the whole run: building and compiling the instrument, and McStas.
        """
        return None if self.finished is None else self.finished - self.started

    def finish(self) -> "Run":
        """Record that the run has finished now, and save the record into the run folder.

        Returns
        -------
        :
            The run with the time now as ``finished``. This run is not changed.
        """
        finished = replace(self, finished=datetime.now().astimezone())
        finished.save()
        return finished

    def save(self) -> Path:
        """Write the record into the run folder, which must exist.

        Returns
        -------
        :
            The file written, ``folder/run.json``.
        """
        record = asdict(self) | {
            "run_type": self.run_type.value,
            "started": self.started.isoformat(),
            "finished": None if self.finished is None else self.finished.isoformat(),
        }
        del record["folder"]  # where the record is found, so a run folder can be moved
        path = self.folder / RUN_FILE
        path.write_text(json.dumps(record, indent=2) + "\n")
        return path

    @classmethod
    def load(cls, folder: str | Path) -> "Run":
        """Read the record of a run folder.

        Parameters
        ----------
        folder:
            The run folder.

        Returns
        -------
        :
            The run, with ``folder`` as its folder.

        Raises
        ------
        FileNotFoundError
            If the folder has no ``run.json``: it is not a finished run.
        """
        folder = Path(folder)
        record = json.loads((folder / RUN_FILE).read_text())
        finished = record.get("finished")  # not in a record from before it was recorded
        record |= {
            "run_type": RunType(record["run_type"]),
            "started": datetime.fromisoformat(record["started"]),
            "finished": None if finished is None else datetime.fromisoformat(finished),
        }
        return cls(folder=folder, **record)


def recorded_instr_params(
    defaults: Iterable[InstrParam], given: Mapping[str, object], variable: str | None = None
) -> dict[str, float | str]:
    """Return the value of every instrument parameter in a run, as its record keeps them.

    Parameters
    ----------
    defaults:
        The instrument parameters of the instrument, with the values it has.
    given:
        The instrument parameters given to the run, by McStas name.
    variable:
        The instrument parameter of a scan, which is left out: its values are recorded apart.

    Returns
    -------
    :
        The parameters given, and the rest as the instrument has them: a number as a float,
        also a NumPy one, and anything else as its string.

    Examples
    --------
    >>> defaults = [InstrParam("omega", 0), InstrParam("b_rot", -0.4), InstrParam("mode", "FF")]
    >>> recorded_instr_params(defaults, {"omega": 30})
    {'omega': 30.0, 'b_rot': -0.4, 'mode': 'FF'}
    >>> recorded_instr_params(defaults, {}, variable="omega")
    {'b_rot': -0.4, 'mode': 'FF'}
    """
    used = {param.name: param.value for param in defaults} | dict(given)
    return {
        name: float(value) if isinstance(value, Real) else str(value)
        for name, value in used.items()
        if name != variable
    }


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


def new_run_folder(runs_folder: Path, folder_name: str | None, label: str) -> Path:
    """Return the folder a run writes into, and clear the way for it.

    The run folder itself is not created, since McStas refuses to write into an existing folder.
    ``runs_folder`` is created with its parents, since McStas does not create them.

    Parameters
    ----------
    runs_folder:
        The folder of the runs.
    folder_name:
        Name of the run folder. A folder of that name from an earlier run is removed, so the run
        replaces it. If ``None``, a new name from the time now and ``label`` is used, so each
        run keeps its own output.
    label:
        What the run does, for the new name, e.g. ``"run"``, ``"scan_omega"``,
        ``"run_to_sample"``.

    Returns
    -------
    :
        ``runs_folder / folder_name``.

    Raises
    ------
    ValueError
        If ``folder_name`` is not a plain folder name of letters, digits, ``.``, ``_`` and ``-``.
    """
    if folder_name is None:
        folder_name = _new_name(runs_folder, label)
    _check_name(folder_name)
    path = runs_folder / folder_name
    shutil.rmtree(path, ignore_errors=True)
    runs_folder.mkdir(parents=True, exist_ok=True)
    return path


def list_runs(runs_folder: str | Path) -> list[Run]:
    """List the finished runs, newest first.

    A run is finished when its folder has a ``run.json``. Any other folder is left out: an
    unfinished run, or something else in ``runs_folder``.

    Parameters
    ----------
    runs_folder:
        The folder of the runs.

    Returns
    -------
    :
        The runs, by the time they started, the newest first. Empty if ``runs_folder`` does not
        exist.
    """
    runs_folder = Path(runs_folder)
    if not runs_folder.is_dir():
        return []
    runs = [Run.load(d) for d in runs_folder.iterdir() if (d / RUN_FILE).is_file()]
    return sorted(runs, key=lambda run: run.started, reverse=True)


def latest_run(
    runs_folder: str | Path, run_type: RunType | None = None, instrument: str | None = None
) -> Run:
    """Return the newest finished run of a type and of an instrument.

    Parameters
    ----------
    runs_folder:
        The folder of the runs.
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
    for run in list_runs(runs_folder):
        if run_type in (None, run.run_type) and instrument in (None, run.instrument):
            return run
    wanted = [
        f"{what} {value}"
        for what, value in (("type", run_type), ("instrument", instrument))
        if value is not None
    ]
    of = f" of {' and '.join(wanted)}" if wanted else ""
    raise FileNotFoundError(f"no finished run{of} in {runs_folder}")


def monitor_file(run: Run) -> Path:
    """Return the file with the incident monitor that normalizes a run.

    A run from the sample has no monitor in its own file, so it is normalized by the run to the
    sample whose beam it started from, its ``beam_run``.

    Parameters
    ----------
    run:
        The run.

    Returns
    -------
    :
        For a run from the sample, the ``mccode.h5`` of its ``beam_run``, in the same runs
        folder. Otherwise its own ``mccode.h5``.

    Raises
    ------
    FileNotFoundError
        If a run from the sample has no ``beam_run``, or that run is not there any more.
    """
    if run.run_type not in (RunType.RUN_FROM_SAMPLE, RunType.SCAN_FROM_SAMPLE):
        return run.data_file
    if run.beam_run is None:
        raise FileNotFoundError(f"{run.folder}: no record of the run to the sample it started from")
    return Run.load(run.folder.parent / run.beam_run).data_file


def clean_runs(
    runs_folder: str | Path, keep: int = 5, run_type: RunType | None = None
) -> list[Path]:
    """Remove all but the newest finished runs.

    Only a finished run is ever removed: a folder without a ``run.json``, an unfinished run or
    anything else in ``runs_folder``, is left as it is.

    Parameters
    ----------
    runs_folder:
        The folder of the runs.
    keep:
        Number of newest runs to keep.
    run_type:
        The type of the runs to clean: the newest ``keep`` of that type are kept, its older ones
        removed, and the runs of every other type left as they are. If ``None``, the runs of
        all types are cleaned together, but a ``run_to_sample`` is never removed: its beam
        file is the expensive simulation every run from the sample starts from, and its
        ``mccode.h5`` has the incident monitor that normalizes them. It still counts toward
        ``keep``, so fewer other runs may be kept. Clean those with
        ``run_type=RunType.RUN_TO_SAMPLE``.

    Returns
    -------
    :
        The removed folders, the newest first.
    """
    runs = list_runs(runs_folder)
    if run_type is None:
        old = [run for run in runs[keep:] if run.run_type is not RunType.RUN_TO_SAMPLE]
    else:
        old = [run for run in runs if run.run_type is run_type][keep:]
    for run in old:
        shutil.rmtree(run.folder)
    return [run.folder for run in old]
