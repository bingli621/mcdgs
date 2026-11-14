from datetime import UTC, datetime, timedelta

import pytest

from mcdgs.mcstas.runs import Run, RunType

STARTED = datetime(2026, 10, 22, 12, tzinfo=UTC)


@pytest.fixture
def finished_run():
    """finished_run(folder, minute, run_type=RunType.RUN, **fields) -> the `Run` of `folder`,
    made a fake finished run: the folder with its `run.json`, started `minute` minutes after
    `STARTED`, of the instrument "TRex" unless `fields` names another."""

    def make(folder, minute, run_type=RunType.RUN, **fields):
        folder.mkdir(parents=True)
        fields = {"instrument": "TRex", "instr_params": {}} | fields
        run = Run(folder, run_type=run_type, started=STARTED + timedelta(minutes=minute), **fields)
        run.save()
        return run

    return make
