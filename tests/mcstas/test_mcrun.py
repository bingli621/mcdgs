import os
import stat
from pathlib import Path

import pytest

from mcdgs.instruments.trex.model import TRex
from mcdgs.mcstas import build_mcstas
from mcdgs.mcstas.mcrun import LOG_FILE, instrument_parameters, mcrun_command, run_mcrun


def test_the_command_has_the_settings_then_the_instrument_then_its_parameters():
    command = mcrun_command(
        "trex.instr",
        Path("runs/a"),
        {"omega": 30.0, "run_to_mcpl": '"beam.mcpl"'},
        ncount=1e5,
        seed=3,
        mpi="auto",
        custom_flags=["--gravitation"],
    )
    assert command == [
        "mcrun",
        "-c",
        "--format=NeXus",
        "-n",
        "100000",
        "--mpi=auto",
        "--seed=3",
        "-d",
        "runs/a",
        "--gravitation",
        "trex.instr",
        "omega=30.0",
        "run_to_mcpl=beam.mcpl",  # without the quotes of the C string, as the shell gives it
    ]


def test_without_ncount_seed_and_mpi_mcrun_keeps_its_defaults():
    command = mcrun_command("trex.instr", Path("a"), {}, nexus=False)
    assert command == ["mcrun", "-c", "-d", "a", "trex.instr"]


def test_the_parameters_are_those_mcstasscript_would_give_mcrun(tmp_path):
    inst = build_mcstas(TRex(), input_path=tmp_path)
    inst.set_parameters(omega=30)
    inst.run_to("sample_pos", filename='"beam.mcpl"')
    parameters = instrument_parameters(inst)
    assert list(parameters) == [p.name for p in inst.parameters]
    assert parameters["omega"] == 30 and parameters["run_to_mcpl"] == '"beam.mcpl"'


@pytest.fixture
def fake_mcrun(tmp_path, monkeypatch):
    """A `mcrun` on the PATH that makes the folder after `-d`, prints its arguments and exits
    with the code in $FAKE_EXIT."""
    bin_folder = tmp_path / "bin"
    bin_folder.mkdir()
    script = bin_folder / "mcrun"
    script.write_text(
        "#!/bin/sh\n"
        'while [ "$1" != "-d" ]; do shift; done\n'
        'mkdir "$2"\n'
        'echo "ran in $(pwd)"\n'
        'exit "${FAKE_EXIT:-0}"\n'
    )
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", f"{bin_folder}{os.pathsep}{os.environ['PATH']}")


@pytest.mark.usefixtures("fake_mcrun")
def test_mcrun_runs_in_the_build_folder_and_its_output_is_kept_in_the_run_folder(tmp_path):
    folder = tmp_path / "run"
    run_mcrun(mcrun_command("x.instr", folder, {}), cwd=tmp_path, folder=folder)
    assert (folder / LOG_FILE).read_text().strip() == f"ran in {tmp_path.resolve()}"


@pytest.mark.usefixtures("fake_mcrun")
def test_a_failed_mcrun_raises_with_its_output_and_keeps_the_log(tmp_path, monkeypatch):
    monkeypatch.setenv("FAKE_EXIT", "3")
    folder = tmp_path / "run"
    with pytest.raises(RuntimeError, match=r"exit code 3, see .*mcrun.log:\nran in"):
        run_mcrun(mcrun_command("x.instr", folder, {}), cwd=tmp_path, folder=folder)
    assert (folder / LOG_FILE).is_file()


def test_without_mcrun_on_the_path_it_says_where_to_run(tmp_path):
    with pytest.raises(FileNotFoundError, match="pixi"):
        run_mcrun(["no-such-mcrun"], cwd=tmp_path, folder=tmp_path / "run")
