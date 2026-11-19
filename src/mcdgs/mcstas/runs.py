"""The runs of a :class:`.Simulation`: their records.

Each run writes into its own folder in the runs folder (:attr:`.Workspace.runs_folder`).

A :class:`Run` says what one run did: which instrument, with which :class:`RunType` and which
values of the instrument parameters, and when it started and finished. It is saved as
``run.json`` in the run folder once McStas has finished (:meth:`Run.finish`), so a folder with
a ``run.json`` is a finished run, and it is read back from there, also in a later session and
without the instrument.

- :func:`recorded_instr_params` gives the instrument parameters of a run as its record keeps
  them.
- :func:`monitor_file` gives the file a run is normalized by.

The :class:`.Workspace` picks the folder of a new run, finds the finished runs and removes older
ones.

Nothing here knows McStas, an :class:`.Instrument` or a :class:`.Workspace`: what a run did is
read from its record, never from the name of its folder.
"""

import json
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timedelta
from enum import StrEnum
from numbers import Real
from pathlib import Path

from mcdgs.core import InstrParam

RUN_FILE = "run.json"  # the record of a run, in its folder
DATA_FILE = "mccode.h5"  # the McStas NeXus file of a run
BEAM_FILE = "beam.mcpl.gz"  # the beam a run to the sample dumps


class RunType(StrEnum):
    """What a run does: the method of :class:`.Simulation` that made it.

    ``DIAGNOSE`` is a run of :meth:`.McStasScriptDiagnostics.run`, ``ANALYSIS`` one of
    :meth:`.McStasScriptDiagnostics.run_analysis`.
    """

    RUN = "run"
    DIAGNOSE = "diagnose"
    ANALYSIS = "analysis"
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
