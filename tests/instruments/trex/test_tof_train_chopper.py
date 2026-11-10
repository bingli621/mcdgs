"""Standalone McStas checks for the train-aware TREX disk chopper."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

COMPONENT_DIR = Path(__file__).parents[3] / "src/mcdgs/instruments/trex/tof_train_components"

pytestmark = pytest.mark.skipif(shutil.which("mcrun") is None, reason="McStas is not on PATH")


INSTRUMENT = r"""DEFINE INSTRUMENT train_chopper(
    int N_trains_par=7,
    double nu=0,
    double delay=0,
    double park_angle=0,
    double slit_lo=0,
    double slit_hi=359.9)

DECLARE
%{
double *tof_train_offsets;
double *tof_train_weights;
long tof_train_arrived;
long tof_train_sent;
long tof_train_rays;
int N_trains;
%}

USERVARS
%{
double *t_offset;
double *p_trains;
double p_last_time_manipulation;
int adaptive_N;
int N_active;
long *train_total_arrived;
long *train_total_N_sent;
long *train_total_rays_sent;
%}

INITIALIZE
%{
N_trains = N_trains_par;
if (N_trains < 1) exit(EXIT_FAILURE);
tof_train_offsets = malloc(N_trains * sizeof(double));
tof_train_weights = malloc(N_trains * sizeof(double));
if (!tof_train_offsets || !tof_train_weights) exit(EXIT_FAILURE);
%}

TRACE
COMPONENT train_init = Arm()
AT (0, 0, 0) ABSOLUTE
EXTEND %{
t_offset = tof_train_offsets;
p_trains = tof_train_weights;
p_last_time_manipulation = 0;
adaptive_N = N_trains;
N_active = 0;
train_total_arrived = &tof_train_arrived;
train_total_N_sent = &tof_train_sent;
train_total_rays_sent = &tof_train_rays;
%}

COMPONENT source = ESS_butterfly(
    sector="W", beamline=7, Lmin=2, Lmax=2.1,
    focus_xw=0.08, focus_yh=0.08, dist=2, n_pulses=1,
    target_tsplit=0)
AT (0, 0, 0.1) ABSOLUTE
COMPONENT chopper = NXdisk_chopper(
    slit_edges={slit_lo, slit_hi}, n_edges=2, radius=0.35,
    nu=nu, delay=delay, park_angle=park_angle, abs_out=1)
AT (0, 0, 1) ABSOLUTE
COMPONENT detector = PSD_monitor(
    nx=1, ny=1, filename="detector.dat", xwidth=0.08, yheight=0.08,
    restore_neutron=1)
AT (0, 0, 1) ABSOLUTE

FINALLY
%{
free(tof_train_offsets);
free(tof_train_weights);
%}
END
"""


def _run(folder: Path, stock: bool = False, **parameters: float) -> float:
    folder.mkdir(parents=True)
    instrument = folder / "train_chopper.instr"
    text = INSTRUMENT
    if stock:
        text = text.replace(
            "slit_edges={slit_lo, slit_hi}, n_edges=2, radius=0.35,\n"
            "    nu=nu, delay=delay, park_angle=park_angle, abs_out=1",
            "theta_0=slit_hi-slit_lo, nslit=1, radius=0.35, yheight=0.35,\n"
            "    nu=nu, delay=delay, abs_out=1",
        ).replace("NXdisk_chopper", "DiskChopper")
    instrument.write_text(text)
    command = [
        "mcrun",
        "-n",
        "4000",
        "--seed=12345",
        "-d",
        str(folder / "run"),
        "-I",
        str(COMPONENT_DIR),
        str(instrument),
        *[f"{name}={value}" for name, value in parameters.items()],
    ]
    result = subprocess.run(command, cwd=folder, text=True, capture_output=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    data = folder / "run" / "detector.dat"
    for line in data.read_text().splitlines():
        if line.startswith("# values:"):
            return float(line.split()[2])
    raise AssertionError(data.read_text())


def test_full_open_and_full_block_match_stock_chopper(tmp_path):
    train_open = _run(tmp_path / "train-open", slit_lo=0, slit_hi=359.9, nu=14)
    stock_open = _run(tmp_path / "stock-open", stock=True, slit_lo=0, slit_hi=359.9, nu=14)
    assert train_open > 0
    assert train_open == pytest.approx(stock_open, rel=0.05)

    train_block = _run(tmp_path / "train-block", slit_lo=0, slit_hi=1, park_angle=180)
    assert train_block == pytest.approx(0, abs=1e-12)


def test_rotating_train_chopper_gates_part_of_each_train(tmp_path):
    all_open = _run(tmp_path / "all-open", slit_lo=0, slit_hi=359.9, nu=14)
    partial = _run(tmp_path / "partial", slit_lo=0, slit_hi=20, nu=14)
    stock = _run(tmp_path / "stock", stock=True, slit_lo=0, slit_hi=20, nu=14)
    assert 0.05 < partial / all_open < 0.25
    assert stock / all_open > 0.7


@pytest.mark.parametrize("nu", [14, -14])
def test_signed_rotation_is_supported(tmp_path, nu):
    all_open = _run(tmp_path / f"open-{nu}", slit_lo=0, slit_hi=359.9, nu=nu)
    partial = _run(tmp_path / f"partial-{nu}", slit_lo=0, slit_hi=20, nu=nu)
    assert 0.05 < partial / all_open < 0.85


def test_parked_chopper_uses_park_angle(tmp_path):
    open_ = _run(tmp_path / "park-open", slit_lo=0, slit_hi=20, park_angle=0)
    blocked = _run(tmp_path / "park-blocked", slit_lo=0, slit_hi=20, park_angle=180)
    assert open_ > 0
    assert blocked == pytest.approx(0, abs=1e-12)
