"""The record of a finished run, kept in its run folder.

A :class:`Run` says what one run of a :class:`.Simulation` did: which instrument, in which
:class:`RunType`, with which values of the instrument parameters. It is saved as ``run.json`` in
the run folder once McStas has finished, so a folder with a ``run.json`` is a finished run, and
it is read back from there, also in a later session and without the instrument.

Nothing here knows McStas, an instrument or a :class:`.Workspace`: a run is found by its folder.
"""

import json
from dataclasses import dataclass
from datetime import datetime
from enum import StrEnum
from pathlib import Path

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

    def save(self) -> Path:
        """Write the record into the run folder, which must exist.

        Returns
        -------
        :
            The file written, ``folder/run.json``.
        """
        record = {
            "instrument": self.instrument,
            "run_type": self.run_type.value,
            "started": self.started.isoformat(),
            "instr_params": self.instr_params,
            "ncount": self.ncount,
            "seed": self.seed,
            "variable": self.variable,
            "values": self.values,
            "beam_run": self.beam_run,
        }
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
        return cls(
            folder=folder,
            instrument=record["instrument"],
            run_type=RunType(record["run_type"]),
            started=datetime.fromisoformat(record["started"]),
            instr_params=record["instr_params"],
            ncount=record["ncount"],
            seed=record["seed"],
            variable=record["variable"],
            values=record["values"],
            beam_run=record["beam_run"],
        )
