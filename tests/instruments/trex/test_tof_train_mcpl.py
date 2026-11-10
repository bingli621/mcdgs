"""Standalone train-to-MCPL checks using the installed McStas runtime."""

import shutil
import subprocess

import pytest
from test_tof_train_chopper import COMPONENT_DIR

mcpl = pytest.importorskip("mcpl")
pytestmark = pytest.mark.skipif(shutil.which("mcrun") is None, reason="McStas is not on PATH")


def run(folder, text, count, train=True, mpi=False):
    instrument = folder / "test.instr"
    instrument.write_text(text)
    result = subprocess.run(
        [
            "mcrun",
            "-y",
            "-n",
            str(count),
            "--seed=123",
            "-d",
            str(folder / "run"),
            "--mpi=2" if mpi else "--no-mpi",
            *(["-I", str(COMPONENT_DIR)] if train else []),
            str(instrument),
        ],
        cwd=folder,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def particles(path):
    with mcpl.MCPLFile(str(path)) as file:
        return [(p.time, p.weight, p.x, p.polx) for p in file.particles]


@pytest.mark.parametrize("weight_mode", [1, 2])
def test_train_output_and_mcpl_input_roundtrip(tmp_path, weight_mode):
    run(
        tmp_path,
        r"""DEFINE INSTRUMENT test()
DECLARE
%{
double offsets[3]; double weights[3];
%}
USERVARS
%{
double *t_offset; double *p_trains; double p_last_time_manipulation; int N_active;
%}
TRACE
COMPONENT source = Arm() AT (0,0,0) ABSOLUTE
EXTEND %{
  x=0.12; vx=0; vy=0; vz=1000; t=0.01; p=8; sx=0.5;
  offsets[0]=0.001; offsets[1]=0.004; offsets[2]=0.009;
  weights[0]=1; weights[1]=3; weights[2]=999;
  t_offset=offsets; p_trains=weights; p_last_time_manipulation=4; N_active=2;
%}
COMPONENT first = MCPL_output(filename="first", doubleprec=1, polarisationuse=1, weight_mode=MODE)
AT (0,0,0) ABSOLUTE
COMPONENT second = MCPL_output(filename="second", doubleprec=1, polarisationuse=1, weight_mode=MODE)
AT (0,0,0) ABSOLUTE
END
""".replace("MODE", str(weight_mode)),
        2,
    )
    first = particles(tmp_path / "run/first.mcpl.gz")
    second = particles(tmp_path / "run/second.mcpl.gz")
    expected_weight = 2 if weight_mode == 2 else 4
    assert len(first) == 4
    assert first == second
    for i in (0, 2):
        assert first[i] == pytest.approx((11, expected_weight, 12, 0.5))
        assert first[i + 1] == pytest.approx((14, 3 * expected_weight, 12, 0.5))

    # The standard reader must consume the files as ordinary neutron events.
    replay = tmp_path / "replay"
    replay.mkdir()
    run(
        replay,
        r"""DEFINE INSTRUMENT test()
TRACE
COMPONENT source = MCPL_input(filename="INPUT", polarisationuse=1)
AT (0,0,0) ABSOLUTE
COMPONENT writer = MCPL_output(filename="replayed", doubleprec=1, polarisationuse=1, weight_mode=2)
AT (0,0,0) ABSOLUTE
END
""".replace("INPUT", str(tmp_path / "run/first.mcpl.gz")),
        4,
        train=False,
    )
    actual = particles(replay / "run/replayed.mcpl.gz")
    # Input normalizes mode 1 by the stored initial ray count.
    divisor = 2 if weight_mode == 1 else 1
    assert len(actual) == len(first)
    for got, original in zip(actual, first):
        assert got == pytest.approx((original[0], original[1] / divisor, original[2], original[3]))


@pytest.mark.skipif(
    not shutil.which("mpirun") or not shutil.which("mpicc"), reason="MPI not available"
)
def test_mpi_train_mcpl_normalization_matches_serial(tmp_path):
    source = r"""DEFINE INSTRUMENT test()
DECLARE
%{
double offsets[2], weights[2], source_normalization;
%}
USERVARS
%{
double *t_offset;
double *p_trains;
double p_last_time_manipulation;
int N_active;
%}
INITIALIZE
%{
source_normalization = 1.0 / mcget_ncount();
%}
TRACE
COMPONENT source = Arm() AT (0,0,0) ABSOLUTE
EXTEND %{
  vz=1000; p=8.0*source_normalization; t=0.01;
  offsets[0]=0.001; offsets[1]=0.004; weights[0]=1; weights[1]=3;
  t_offset=offsets; p_trains=weights; p_last_time_manipulation=4; N_active=2;
%}
COMPONENT writer = MCPL_output(filename="events", doubleprec=1)
AT (0,0,0) ABSOLUTE
END
"""
    replay = r"""DEFINE INSTRUMENT test()
TRACE
COMPONENT source = MCPL_input(filename="INPUT") AT (0,0,0) ABSOLUTE
COMPONENT writer = MCPL_output(filename="replayed", doubleprec=1, weight_mode=2)
AT (0,0,0) ABSOLUTE
END
"""
    totals = []
    for mpi in (False, True):
        folder = tmp_path / ("mpi" if mpi else "serial")
        folder.mkdir()
        run(folder, source, 100, mpi=mpi)
        path = folder / "run/events.mcpl.gz"
        with mcpl.MCPLFile(str(path)) as file:
            assert file.nparticles == 200
            assert file.stat_sum["initial_ray_count"] == pytest.approx(100)
        replay_dir = folder / "replay"
        replay_dir.mkdir()
        run(replay_dir, replay.replace("INPUT", str(path)), 200, train=False)
        events = particles(replay_dir / "run/replayed.mcpl.gz")
        assert len(events) == 200
        totals.append(sum(event[1] for event in events))
    assert totals == pytest.approx([8, 8])


def test_long_command_line_header_is_safely_truncated(tmp_path):
    parameters = ", ".join(f"double parameter_with_long_name_{i}=1" for i in range(120))
    instrument = r"""DEFINE INSTRUMENT test(PARAMETERS)
USERVARS
%{
double *t_offset;
double *p_trains;
double p_last_time_manipulation;
int N_active;
%}
TRACE
COMPONENT source = Arm() AT (0,0,0) ABSOLUTE
EXTEND %{
  N_active=0; p_last_time_manipulation=0;
%}
COMPONENT writer = MCPL_output(filename="events") AT (0,0,0) ABSOLUTE
END
""".replace("PARAMETERS", parameters)
    run(tmp_path, instrument, 1)
    with mcpl.MCPLFile(str(tmp_path / "run/events.mcpl.gz")) as file:
        command = file.blobs["mccode_cmd_line"]
        assert len(command) == 2047
        assert b"parameter_with_long_name_0=1" in command
