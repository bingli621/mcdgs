# Experimental component-only TOF trains

The train algorithm derives from McStas commit `a247ebb2f`, prototype
`ODIN_TOF_train7`. The source and display geometry are ported onto McStas
main commit `63a9b7a6f5ecd668094c05af5cef135e5b809d71`, including its TBL
geometry and logo option. The monitor and its mesh/monitor libraries are also
ported from that main commit. The four files `trex-tof-train-monitor_nd-lib.{h,c}`
and `trex-tof-train-geometry-mesh-lib.{h,c}` have train-specific include names
and C identifiers to avoid shadowing or colliding with McStas libraries used by
other components. Their upstream code is otherwise unchanged. They use the
installed McStas runtime and `read_table-lib`; keep all four alongside
`Monitor_nD.comp` on the component search path. The two local
`ESS_butterfly-lib` files and `trex-tof-train-butterfly-geometry.c` must be
resolved together with the local component. The local library retains train7's
factorized brightness functions, so time-independent work is done once per ray;
the uniquely named geometry helper uses the current 28-argument interface.
Put this directory first in the McStas component search path; do not mix this
source with the stock library.

## Using the TREX experiment

The optimization is disabled by default. To enable it for a simulation:

```python
from mcdgs.instruments.trex.model import TRex
from mcdgs.mcstas.simulation import Simulation

simulation = Simulation(TRex(), tof_train=True)
```

Alternatively, use `build_mcstas(TRex(), tof_train=True)` from
`mcdgs.mcstas.builder` to obtain the McStasScript instrument. Generated files
remain in the normal input directory; the builder configures component
discovery and McStas `SEARCH` for this directory. The default maximum train
size is `N_trains_par=300` and the target is `target_tsplit=10`; they are McStas
instrument parameters. The pre-sample monitor supplies adaptive feedback.

This experiment currently accepts only the reviewed sample-free TREX layout.
OpenMP, OpenACC, custom extensions and focused emission time are unsupported.
`SPLIT` produces a warning: it must be after the last chopper, with no absorbing
components that modify the train in downstream branches; shared train arrays
and time-dependent components may otherwise produce incorrect results.
`Simulation.run_to_sample()` writes the expanded train
to MCPL; downstream runs read ordinary particles without train components.
MPI ranks have independent feedback,
but MPI execution has not been validated.

## Instrument integration contract

The instrument must define `N_trains` as a positive integer expression visible
to the source, e.g. `#define N_trains INSTRUMENT_GETPAR(N_trains_par)`, and
declare the following particle variables in `USERVARS`:

```c
double *t_offset;
double *p_trains;
double p_last_time_manipulation;
int adaptive_N;
int N_active;
long *train_total_arrived;
long *train_total_N_sent;
long *train_total_rays_sent;
```

Before the source runs, initialize the two pointers for every particle to
distinct writable arrays of at least `N_trains` doubles. Their lifetime must
cover the complete particle trace, including any split/restore behavior. Set
the three `train_total_*` pointers to persistent, zero-initialized `long`
counters, and initialize the remaining scalar variables. Do not use one shared array when
particles can run concurrently or coexist. Allocate and release storage with a
strategy compatible with the chosen runtime and its particle copies; a plain
per-ray allocation with no matching cleanup leaks when particles are absorbed.
The source resets the train weight sum for each new ray.
An all-zero train is absorbed at the source; the gate and read macros also
absorb when their stored total is zero rather than dividing by it.

The source initially generates `N_trains` time samples. Its `target_tsplit`
(default 10) controls the adaptive sample count: after feedback it uses
`ceil(target_tsplit * *train_total_N_sent / *train_total_arrived)`, clamped to 1 through
`N_trains`; no arrivals yet means use `N_trains`. To enable feedback, mark
**one** downstream train-aware monitor `adaptive_target=1`. The source adds
each ray's sampled count to `*train_total_N_sent` and increments
`*train_total_rays_sent` before emission, including for rays later absorbed.
The target adds the number of active samples reaching its geometry to
`*train_total_arrived`, before monitor acceptance, in either split or unsplit
mode. Point every particle in a serial worker at the same counters; initialize
them once per run, not once per ray. Use independent counters per MPI process.
These plain C updates are not atomic: do not share them across simultaneous
threads or GPU work items. Per-worker arrays and counters are required, and
parallel execution has not been validated.

`Monitor_nD` defaults to `tsplit=-1`: it inspects parsed coordinates and
splits when a time coordinate or signal is present. Set `tsplit=1` explicitly
for user variables or other expressions indirectly depending on time, and
`tsplit=0` for a time-independent monitor using the aggregate ray weight.
In split mode a monitor with `absorb` removes only the rejected members and
updates the remaining total weight, including with `restore_neutron=1`.
Every time-dependent component between source and monitor must also understand
the train representation; an ordinary chopper only sees the base time and is
not correct. The accompanying integration must supply a train-aware chopper.

The local `MCPL_output` is a train-aware save component adapted from
`mcstas-comps/misc/MCPL_output.comp` at McCode main commit
`63a9b7a6f5ecd668094c05af5cef135e5b809d71`. It retains the standard
`filename`, `weight_mode`, `verbose`, `polarisationuse`, `doubleprec`,
`userflag`, `userflagcomment`, and `buffermax` parameters and uses the
installed modern MCPL library and runtime. For each active member it writes
an ordinary neutron with time `t + t_offset[i]` and weight
`p * p_trains[i] / p_last_time_manipulation`, followed by the standard
MCPL weight-mode scaling. It restores the incoming time and weight and does
not implicitly serialize train state; the ordinary explicit `userflag` option
still stores the user-selected 32-bit value. An empty or zero-total train
produces no entries. The output can be read by standard `MCPL_input`, but
train correlation and the train representation cannot be reconstructed from
the resulting ordinary particle records. OpenACC is explicitly rejected. A standalone
two-rank MPI regression checks merged `initial_ray_count` and replay intensity
against serial for the source's initialization-time normalization; full TREX
MPI execution remains unvalidated. The embedded command-line metadata is
truncated to 2047 bytes when necessary.

## Limitations

* A split monitor calls the standard Monitor_nD once per time sample. Its `N`
  and sum-of-squared-weights error estimate therefore count correlated samples
  as independent. Intensities are weighted by `1/N`, but the reported error is
  **not** the uncertainty across independent source rays. Use unsplit monitors
  for independent-ray errors or implement a per-original-ray, per-bin
  accumulation before relying on uncertainty estimates. This prototype does
  not provide that accumulation.
* Time focusing (`tfocus_width > 0`) is rejected at initialization: the
  original prototype's whole-train rejection was biased.
* GPU/OpenACC, MPI, concurrent particles, neutron splitting, and restoring
  neutron state have not been validated with this component-only integration.
