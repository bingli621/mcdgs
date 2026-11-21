"""Run McStas with mcrun, on the ``.instr`` that McStasScript writes.

- :func:`instrument_parameters` gives the value of every instrument parameter of a McStasScript
  instrument, as mcrun takes them.
- :func:`mcrun_command` builds the mcrun command of one run.
- :func:`run_mcrun` runs it, and keeps mcrun's output in ``mcrun.log`` in the run folder.

Nothing here knows an :class:`.Instrument`, a :class:`.Simulation` or a :class:`.Run`.
"""

import shutil
import subprocess
from collections.abc import Mapping, Sequence
from pathlib import Path

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


def _argument(name: str, value: object) -> str:
    """Return ``name=value`` as mcrun takes it.

    A string value is a C string, in quotes, to McStasScript; mcrun takes it without them, as the
    shell gives it when McStasScript runs mcrun.
    """
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
    nexus: bool = True,
    custom_flags: Sequence[str] = (),
    mcrun: str = "mcrun",
) -> list[str]:
    """Return the mcrun command of one run, which compiles the instrument and runs it.

    The arguments come in the order McStasScript gives them.

    Parameters
    ----------
    instr_file:
        The ``.instr`` file, relative to the folder mcrun runs in.
    folder:
        The run folder mcrun writes into, which must not exist yet.
    parameters:
        The value of every instrument parameter, by McStas name, as
        :func:`instrument_parameters` gives them.
    ncount:
        Number of neutrons. If ``None``, mcrun's default, as from a beam file, which has its own
        number of neutrons.
    seed:
        Random seed. If ``None``, McStas picks one.
    mpi:
        Number of MPI processes, or ``"auto"``. If ``None``, no MPI.
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
