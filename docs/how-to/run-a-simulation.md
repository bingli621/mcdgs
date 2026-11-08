# Run a McStas simulation

A [`Simulation`](../reference/mcstas.md) holds an instrument and the folder its runs write to. The
examples use `TRex`, as in [Install and first run](../getting-started.md). For the timing of the
choppers alone, [tof](simulate-the-choppers-with-tof.md) is much faster.

```python
import numpy as np

from mcdgs.instruments.trex.model import TRex
from mcdgs.mcstas.simulation import Simulation

simulation = Simulation(TRex())
```

## One run

Keywords that are not settings of the run are McStas instrument parameters. Those you leave out
keep the instrument's values.

```python
h5_file = simulation.run(ncount=1e6, seed=1, bender_in=1)
```

`run` uses all cores through MPI by default (`mpi="auto"`).

## A scan

`scan` runs once per value of an instrument parameter and merges the results into one file, with
one entry per value. This is the [bender scan](../instruments/trex/bender-scan.md#the-scan), with fewer
neutrons and at the sample position:

```python
h5_file = simulation.scan(
    "b_rot",
    np.linspace(-0.8, 0.1, 19),
    ncount=1e5,  # per value
    seed=1,
    bender_in=1,
)
```

The instrument is compiled once, then the values run in parallel as single-process McStas jobs.
`n_jobs` limits how many run at a time. With a `seed`, value `n` uses `seed + n + 1`, so the values
do not repeat each other's random numbers. A name that is not an instrument parameter raises a
`KeyError` before anything runs.

## Find and clean up runs

Each run writes a new folder in `mcstas_output/`, named after its start time and what it did, such
as `20261022-143005_scan_b_rot`.

```python
simulation.list_runs()  # folder names, newest first
simulation.latest_run("*_scan_b_rot")  # the newest folder matching a pattern
simulation.clean_runs(keep=5)  # remove all but the 5 newest
```

`clean_runs` never removes a `run_to_sample` run, nor a folder without a `mccode.h5`.

!!! warning "A `folder_name` is replaced on every run"

    Giving a run a name, as in `simulation.run("bender_out", ...)`, makes the output easy to find
    from a plotting script. An earlier run of that name is removed first.

## A split run at the sample

Simulating the primary spectrometer is the expensive part. A split run does it once, dumps the beam
at the component `sample_pos` to an MCPL file, and starts every later run from there. The
instrument needs a component named `sample_pos`, and a detector after it for the second part, as
T-Rex has.

```python
from mcdgs.instruments.trex.model import TRex
from mcdgs.samples import aluminum

simulation = Simulation(TRex(sample=aluminum()))
beam_file = simulation.run_to_sample(ncount=1e8)  # once, with many neutrons

simulation.run_from_sample()  # part 2, also from a later session
simulation.scan_from_sample("omega", [0, 10, 20])  # part 2, per goniometer angle
```

The beam file is found again by a later `Simulation` on the same `output_dir`. The number of
neutrons in part 2 is the number in the file.

- Make a new beam file after changing a component up to `sample_pos`, or an instrument parameter
  they use.
- Part 2 runs in the frame of `sample_pos`, so anything after it must be placed `relative` to it,
  not in the lab frame.

T-Rex places its `beamstop` and `banana_detector` `relative` to `sample_pos`. Its
`beamstop_monitor` comes before `sample_pos`, so it is counted in part 1, not part 2.
