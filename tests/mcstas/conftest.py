import os

import pytest


@pytest.fixture
def finished_run():
    """finished_run(folder, mtime, beam=False) -> `folder`, made a fake finished run: a
    `mccode.h5` (and a `beam.mcpl.gz` with `beam`), the `mccode.h5` last changed at `mtime`."""

    def make(folder, mtime, beam=False):
        folder.mkdir(parents=True)
        (folder / "mccode.h5").touch()
        if beam:
            (folder / "beam.mcpl.gz").touch()
        os.utime(folder / "mccode.h5", (mtime, mtime))
        return folder

    return make
