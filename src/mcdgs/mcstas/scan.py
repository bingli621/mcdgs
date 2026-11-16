"""Run a McStas instrument once per point of a scan, in parallel, into one scan file.

:func:`run_points_in_parallel` compiles the instrument once (:func:`compile_instrument`), runs
a copy of it per point as a single-process mcrun job, and merges the outputs into one McStas
NeXus file with an entry per point.

It works on a McStasScript instrument and a folder: nothing here knows an :class:`.Instrument`,
a :class:`.Simulation` or a :class:`.Run`.
"""

import copy
import os
import shlex
import shutil
import subprocess
from collections.abc import Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from mcstasscript.interface.instr import McStas_instr

from mcdgs.mcstas.nexus_io import merge_entries
from mcdgs.mcstas.runs import DATA_FILE


def random_seed() -> int:
    """Return a random seed for mcrun, which does not accept 0 or more than nine digits."""
    return int(np.random.SeedSequence().generate_state(1)[0] % 10**9)


def compile_instrument(inst: McStas_instr) -> None:
    """Write the instrument into its input folder and compile it there, without simulating.

    mcrun compiles (``-c``), then the instrument only prints its info (``--info``), so not even a
    beam file is read. The compile flags are those of ``inst``'s settings (NeXus, OpenACC,
    custom flags), so its runs reuse the executable.

    Raises
    ------
    RuntimeError
        With mcrun's output, if compiling fails.
    """
    inst.write_full_instrument()
    options = inst._run_settings  # McStasScript keeps them private
    mcrun = Path(options.get("executable_path", "")) / options.get("executable", "mcrun")
    command = [str(mcrun), "-c", "--info"]
    if options.get("NeXus"):
        command.append("--format=NeXus")
    if options.get("openacc"):
        command.append("--openacc")
    command += [*shlex.split(options.get("custom_flags", "")), f"{inst.name}.instr"]
    process = subprocess.run(
        command, cwd=inst.input_path, capture_output=True, text=True, check=False
    )
    if process.returncode != 0:
        raise RuntimeError(f"compiling {inst.name} failed:\n{process.stdout}{process.stderr}")


def run_points_in_parallel(
    inst: McStas_instr,
    output_path: Path,
    points: Sequence[Mapping[str, float]],
    seed: int | None = None,
    n_jobs: int | None = None,
) -> Path:
    """Run the instrument once per point and merge the outputs into one scan file.

    The instrument is written and compiled once, first, without simulating. Then every point
    runs a copy of ``inst`` in parallel, with its own values, output folder and random seed.

    Parameters
    ----------
    inst:
        The McStasScript instrument, with its settings.
    output_path:
        Folder to create for the scan.
    points:
        Instrument parameter values (by McStas name) of each point.
    seed:
        Seed of the first point; point ``n`` gets ``seed + n + 1``, so the runs do not repeat
        each other's random numbers. If ``None``, a random seed.
    n_jobs:
        Number of points run at a time. If ``None``, the number of cores.

    Returns
    -------
    :
        ``output_path/mccode.h5``, with one entry per point in the order of ``points``.

    Raises
    ------
    ValueError
        If there are no points, before anything runs.
    RuntimeError
        Naming the point that failed, once the points already running are done. Points not
        started yet are not run, and every point's folder is kept in ``output_path``.
    """
    if not points:
        raise ValueError("a scan needs at least one point")
    compile_instrument(inst)  # once, before any point
    output_path.mkdir()  # holds a folder per point
    if seed is None:
        seed = random_seed()
    folders = [output_path / f"point_{n:04d}" for n in range(len(points))]

    runs = []
    for n, (values, folder) in enumerate(zip(points, folders)):
        run = copy.deepcopy(inst)  # a run's settings and values stay its own
        run.settings(
            output_path=str(folder),
            seed=seed + n + 1,  # mcrun does not accept seed 0
            force_compile=False,  # compiled once, above
        )
        run.set_parameters(**values)
        runs.append(run)

    def run_point(n: int) -> None:
        try:
            runs[n].backengine()  # a separate mcrun process
        except Exception as error:
            raise RuntimeError(f"point {n} {points[n]} failed, see {folders[n]}") from error

    with ThreadPoolExecutor(max_workers=n_jobs or os.cpu_count()) as pool:
        futures = [pool.submit(run_point, n) for n in range(len(runs))]
        for future in as_completed(futures):
            if (error := future.exception()) is not None:
                pool.shutdown(wait=False, cancel_futures=True)  # the exit waits for the running
                raise error

    merged = merge_entries([folder / DATA_FILE for folder in folders], output_path / DATA_FILE)
    for folder in folders:
        shutil.rmtree(folder)
    return merged
