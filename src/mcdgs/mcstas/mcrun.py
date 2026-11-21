"""Run McStas with mcrun, on the ``.instr`` that McStasScript writes.

- :func:`instrument_parameters` gives the value of every instrument parameter of a McStasScript
  instrument, as mcrun takes them.
- :func:`mcrun_command` builds the mcrun command of one run, or of a scan: mcrun's own list
  scan, ``-L``, run in parallel with ``--scan_split``.
- :func:`run_mcrun` runs it, and keeps mcrun's output in ``mcrun.log`` in the run folder.
- :func:`check_scan` checks that a scan file has an entry for every value.
- :func:`random_seed` picks a seed, to record it.

A scan writes one McStas NeXus file with an entry per value, ``entry1``, ``entry2``, ...: mcrun
runs each value in a folder of its own, then merges them and removes the folders. Value ``n``
takes the seed ``seed + n * 1024``.

Nothing here knows an :class:`.Instrument`, a :class:`.Simulation` or a :class:`.Run`.
"""

import shutil
import subprocess
from collections.abc import Mapping, Sequence
from pathlib import Path

import mcstastox
import numpy as np
from mcstasscript.interface.instr import McStas_instr

LOG_FILE = "mcrun.log"  # mcrun's output, in the run folder
LOG_TAIL = 40  # the lines of mcrun's output an error shows


def instrument_parameters(inst: McStas_instr) -> dict[str, object]:
    """Return the value of every instrument parameter of a McStasScript instrument.

    These are the values McStasScript itself would give mcrun: those set with
    ``set_parameters``, and else the instrument's, including the beam file parameters
    ``run_to_mcpl`` and ``run_from_mcpl`` that ``run_to`` and ``run_from`` add.

    Parameters
    ----------
    inst:
        The McStasScript instrument.

    Returns
    -------
    :
        The values by McStas name, in the order of the instrument.

    Raises
    ------
    ValueError
        If a parameter has no value.
    """
    values = {}
    for parameter in inst.parameters:
        if parameter.value is None:
            raise ValueError(f"{inst.name}: the instrument parameter {parameter.name} has no value")
        values[parameter.name] = parameter.value
    return values


def random_seed() -> int:
    """Return a random seed for mcrun, which does not accept 0.

    It is below 10**9, so that the seeds of the values of a scan, ``seed + n * 1024``, stay
    below 2**31 for up to a million values.
    """
    return int(np.random.SeedSequence().generate_state(1)[0] % (10**9 - 1)) + 1


def _argument(name: str, value: object) -> str:
    """Return ``name=value`` as mcrun takes it.

    A string value is a C string, in quotes, to McStasScript; mcrun takes it without them, as the
    shell gives it when McStasScript runs mcrun. A list is the values of a scan, ``-L``, which
    mcrun takes separated by commas.
    """
    if isinstance(value, list):
        return f"{name}={','.join(str(item) for item in value)}"
    if isinstance(value, str) and len(value) >= 2 and value[0] == value[-1] == '"':
        value = value[1:-1]
    return f"{name}={value}"


def mcrun_command(
    instr_file: str,
    folder: Path,
    parameters: Mapping[str, object],
    *,
    ncount: float | None = None,
    seed: int | None = None,
    mpi: int | str | None = None,
    scan_split: int | str | None = None,
    nexus: bool = True,
    custom_flags: Sequence[str] = (),
    mcrun: str = "mcrun",
) -> list[str]:
    """Return the mcrun command of one run or scan, which compiles the instrument and runs it.

    The arguments come in the order McStasScript gives them, and those of a scan after the seed.

    Parameters
    ----------
    instr_file:
        The ``.instr`` file, relative to the folder mcrun runs in.
    folder:
        The run folder mcrun writes into, which must not exist yet.
    parameters:
        The value of every instrument parameter, by McStas name, as
        :func:`instrument_parameters` gives them; of a scan, the values of its parameter as a
        list.
    ncount:
        Number of neutrons. If ``None``, mcrun's default, as from a beam file, which has its own
        number of neutrons.
    seed:
        Random seed. If ``None``, McStas picks one.
    mpi:
        Number of MPI processes, or ``"auto"``. If ``None``, no MPI, as for a scan, whose values
        run in parallel instead.
    scan_split:
        Of a scan, the number of values run at a time, or ``"auto"``, the number of cores less
        one. If ``None``, the run is no scan.
    nexus:
        True writes McStas NeXus, ``mccode.h5``.
    custom_flags:
        Extra mcrun arguments.
    mcrun:
        The mcrun executable.

    Returns
    -------
    :
        The command, one argument per item.

    Examples
    --------
    >>> mcrun_command("small.instr", Path("runs/a"), {"width": 0.1}, ncount=1e4, seed=1)
    ['mcrun', '-c', '--format=NeXus', '-n', '10000', '--seed=1', '-d', 'runs/a', 'small.instr', 'width=0.1']
    >>> mcrun_command("small.instr", Path("a"), {"width": [0.1, 0.2]}, seed=1, scan_split=2)[4:]
    ['-L', '--scan_split=2', '-d', 'a', 'small.instr', 'width=0.1,0.2']
    """
    command = [mcrun, "-c"]  # -c compiles, as McStasScript does on every run
    if nexus:
        command.append("--format=NeXus")
    if ncount is not None:
        command += ["-n", str(int(ncount))]
    if mpi is not None:
        command.append(f"--mpi={mpi}")
    if seed is not None:
        command.append(f"--seed={seed}")
    if scan_split is not None:
        command += ["-L", f"--scan_split={scan_split}"]
    command += ["-d", str(folder), *custom_flags, instr_file]
    return command + [_argument(name, value) for name, value in parameters.items()]


def run_mcrun(command: Sequence[str], cwd: Path, folder: Path) -> None:
    """Run mcrun, and keep its output in ``mcrun.log`` in the run folder.

    Parameters
    ----------
    command:
        The mcrun command, as :func:`mcrun_command` gives it; its first item is found on the
        ``PATH``.
    cwd:
        The folder mcrun runs in, which holds the ``.instr`` file.
    folder:
        The run folder of the command, where the log goes once mcrun has made it.

    Raises
    ------
    FileNotFoundError
        If mcrun is not on the ``PATH``.
    RuntimeError
        With the end of mcrun's output, if mcrun fails.
    """
    executable = shutil.which(command[0])
    if executable is None:
        raise FileNotFoundError(f"{command[0]} is not on the PATH: run in the pixi environment")
    process = subprocess.run(
        [executable, *command[1:]], cwd=cwd, capture_output=True, text=True, check=False
    )
    output = process.stdout + process.stderr
    log = folder / LOG_FILE
    if folder.is_dir():  # mcrun makes it, unless it fails first
        log.write_text(output)
    if process.returncode != 0:
        tail = "\n".join(output.splitlines()[-LOG_TAIL:])
        where = f", see {log}" if folder.is_dir() else ""
        raise RuntimeError(f"mcrun failed with exit code {process.returncode}{where}:\n{tail}")


def check_scan(folder: Path, variable: str, values: Sequence[float]) -> None:
    """Check that the file of a scan has an entry for every value.

    mcrun skips a value whose run fails, or prints no ``Detector:`` line, as an instrument whose
    only monitors are event lists, keeps the folders of the values, and merges the rest. It
    exits without an error even if it skips them all.

    Parameters
    ----------
    folder:
        The run folder of the scan, with its ``mccode.h5``.
    variable, values:
        The instrument parameter of the scan and its values.

    Raises
    ------
    RuntimeError
        If the file has fewer entries than values.
    """
    try:
        with mcstastox.Read(str(folder)) as data:
            n_entries = data.get_number_of_entries()
    except ValueError:  # mcstastox does not open a file without entries
        n_entries = 0
    if n_entries != len(values):
        raise RuntimeError(
            f"{folder}: mcrun ran {n_entries} of the {len(values)} values of {variable}, "
            f"see {folder / LOG_FILE}"
        )
