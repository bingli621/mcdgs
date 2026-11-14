# Run a McStas simulation

A [`Simulation`](../reference/mcstas.md) holds an instrument and the
[workspace](#where-the-files-go) its runs write into. The examples use `TRex`, as in
[Install and first run](../getting-started.md). For the timing of the choppers alone,
[tof](simulate-the-choppers-with-tof.md) is much faster.

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
run = simulation.run(ncount=1e6, seed=1, bender_in=1)
run.data_file  # mccode.h5, in the run folder run.folder
```

`run` uses all cores through MPI by default (`mpi="auto"`). Like every run method, it returns the
[record of the run](#find-and-clean-up-runs), which knows its folder and its McStas NeXus file.

## A scan

`scan` runs once per value of an instrument parameter and merges the results into one file, with
one entry per value. This is the [bender scan](../instruments/trex/bender-scan.md#the-scan), with fewer
neutrons and at the sample position:

```python
scan = simulation.scan(
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

## Where the files go

Everything is written into one folder, the workspace:

| Folder | Contents |
| ------ | -------- |
| `mcdgs_workspace/builds/<instrument name>/` | The generated `.instr` and `.c` files and the compiled instrument. |
| `mcdgs_workspace/runs/<run name>/` | What McStas wrote in one run, and the record of the run. |
| `mcdgs_workspace/results/` | For what is made from the runs, such as figures. |

Without anything else, `mcdgs_workspace/` is in the folder you run Python in. To keep it in one
place, set the environment variable `MCDGS_WORKSPACE` to a folder, or give the simulation a
workspace:

```python
from mcdgs.mcstas.workspace import Workspace

workspace = Workspace.find("~/trex_workspace")
simulation = Simulation(TRex(), workspace)
```

Each run writes a new folder in `runs/`, named after its start time and what it did, such as
`20261022-143005_scan_b_rot`.

!!! warning "A `folder_name` is replaced on every run"

    Giving a run a name, as in `simulation.run("bender_out", ...)`, makes the output easy to find
    from a plotting script. An earlier run of that name is removed first.

## Find and clean up runs

A run that finished leaves its record, `run.json`, in its folder: the instrument, the type of the
run, the number of neutrons, the seed, and the value of every instrument parameter. A folder
without one is not a finished run. The runs are found by that record, also in a later session and
without the instrument:

```python
from mcdgs.mcstas.runs import RunType, clean_runs, latest_run, list_runs

runs_folder = simulation.workspace.runs_folder
list_runs(runs_folder)  # the finished runs, newest first
scan = latest_run(runs_folder, RunType.SCAN)  # the newest scan
scan.folder, scan.data_file  # where it is, and its mccode.h5
scan.variable, scan.values, scan.instr_params  # what it simulated
```

`clean_runs` removes all but the newest runs:

```python
clean_runs(runs_folder, keep=5)  # all but the 5 newest
clean_runs(runs_folder, keep=1, run_type=RunType.RUN_TO_SAMPLE)  # all but the newest beam
```

- Without a `run_type`, it never removes a `run_to_sample` run, whose beam file the
  [split runs](#a-split-run-at-the-sample) start from.
- With a `run_type`, it removes only runs of that type.
- It never removes a folder without a `run.json`.

The build folders only hold what every run makes anew, so `simulation.workspace.clean_builds()`
removes them all.

## A split run at the sample

Simulating the primary spectrometer is the expensive part. A split run does it once, dumps the beam
at the component `sample_pos` to an MCPL file, and starts every later run from there. The
instrument needs a component named `sample_pos`, and a detector after it for the second part, as
T-Rex has.

```python
from mcdgs.instruments.trex.model import TRex
from mcdgs.samples import aluminum

simulation = Simulation(TRex(sample=aluminum()))
to_sample = simulation.run_to_sample("to_sample", ncount=1e8)  # once, with many neutrons
to_sample.beam_file  # beam.mcpl.gz, in its run folder

simulation.run_from_sample()  # part 2, from the beam of the run above
simulation.scan_from_sample("omega", [0, 10, 20])  # part 2, per goniometer angle
```

A later session names the run to the sample whose beam it starts from, which must be in the same
workspace:

```python
simulation = Simulation(TRex(sample=aluminum()))
simulation.scan_from_sample("omega", [0, 10, 20], beam_run="to_sample")
```

The beam of an earlier session is never picked by itself: the newest one in the workspace may be
that of another instrument. The number of neutrons in part 2 is the number in the file.

Part 2 has no incident monitor of its own. Its record names the run to the sample it started from,
and `monitor_file` gives the file of that run:

```python
from mcdgs.mcstas.runs import monitor_file

scan = latest_run(runs_folder, RunType.SCAN_FROM_SAMPLE)
scan.beam_run  # "to_sample"
monitor_file(scan)  # the mccode.h5 of "to_sample"
```

- Make a new beam file after changing a component up to `sample_pos`, or an instrument parameter
  they use.
- Removing a run to the sample leaves the runs from it without their monitor.
- Part 2 runs in the frame of `sample_pos`, so anything after it must be placed `relative` to it,
  not in the lab frame.

T-Rex places its `beamstop` and `banana_detector` `relative` to `sample_pos`. Its
`beamstop_monitor` comes before `sample_pos`, so it is counted in part 1, not part 2.
