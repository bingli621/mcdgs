# Core concepts

## An instrument is a class

An instrument is a subclass of [`Instrument`](reference/core.md) holding two things:

- its **components**, one attribute each, in beam order;
- its **control parameters**, `params`, a frozen dataclass.

The attribute name of a component is its McStas instance name, and the name other components use
to refer to it. `trex.guide_1` is the component `guide_1` in the generated `.instr` file.

## Components

A component is a dataclass named after the McStas component it stands for, such as
`Guide_gravity`. Its fields are the McStas parameters, with McStas's names and defaults. See
[Components](reference/components.md) for the ones that exist.

Every component also has the McStas instance keywords: `at`, `rotated`, `relative`, `when`,
`group`, `jump`, `split` and `extend`. For example, T-Rex's guide 4a only acts with the bender out:

```python
Guide_gravity(..., relative="arm_guide_4", when="bender_in==0")
```

A value can also be the name of a McStas instrument parameter, or an expression of them, as in
McStas. The bender rotates with `rotated=[0, "b_rot", 0]`.

## Four kinds of parameter

| Kind | Where it lives | Changed by | Example |
| ---- | -------------- | ---------- | ------- |
| Design parameter | The fields of a component | Editing the instrument | `guide_1.w1` |
| Control parameter | `instrument.params` | `instrument.update_params(...)` | `rrm=18`, `mode="FF"`, `bender_in=True` |
| McStas instrument parameter | `instrument.instr_params`, computed from the control parameters | A keyword of a run | `nu_mc1`, `delay_mc1`, `bender_in` |
| Simulation parameter | Neither the instrument nor its parameters | A keyword of a run or a build | `simulation.run(ncount=1e6)`, `build_tof(trex, neutrons=1_000_000)` |

The control parameters are the settings of the instrument as it is operated, such as the
repetition rate multiplication `rrm` and the chopper mode `mode` of T-Rex.

The McStas instrument parameters are the ones McStas knows about, and the instrument computes them
from its control parameters. For T-Rex, the frequency `nu_<chopper>` and delay `delay_<chopper>`
of each chopper follow from `lambda0`, `rrm`, `mode` and `t0_offset`. Their values in
`instr_params` are the defaults of a run, and a keyword of a run overrides them for that run only,
without changing the instrument, as in `simulation.run(bender_in=1)`.

The simulation parameters say how to simulate, not what: the number of neutrons, `ncount` of a
McStas run or `neutrons` of a tof source, a random seed, or the number of MPI processes. They
belong to one run, or one tof model, and are never stored with the instrument.

## Coordinates

Positions and rotations are those of McStas: z along the beam, y up, x to the left of the beam.
`at` is in m and `rotated` in deg, about x, then y, then z, both relative to the component named by
`relative`, or to the lab frame if there is none. `instrument.positions()` follows these chains
back to the lab frame.

## From instrument to result

An instrument can be simulated with McStas or with tof. Neither changes the instrument, so the
same object can go to both.

### With McStas

1. [`build_mcstas`](reference/mcstas.md) turns the instrument into a McStasScript instrument. It
   adds the instrument parameters, then the components in beam order.
2. [`Simulation`](reference/mcstas.md) writes it to `mcstas_input/<instrument name>/`, where McStas
   compiles it, and runs it. Every run builds the instrument anew, so a change applies to the next
   run.
3. Each run writes a folder in `mcstas_output/` with a McStas NeXus file, `mccode.h5`.
4. [`read_pixels`](reference/mcstas.md) or `h5py` reads that file. See
   [Read the output](how-to/read-the-output.md).

### With tof

1. [`build_tof`](reference/tof.md) turns the instrument into a tof model: the source, the choppers,
   and a detector at every monitor. It leaves out the optics, and measures each distance in a
   straight line from the source.
2. `model.run()` follows the neutrons in memory and writes no files.
3. `result["sample_monitor"]` holds the neutrons at the monitor at the sample position, with their
   wavelength and time of arrival. See [Simulate the choppers with tof](how-to/simulate-the-choppers-with-tof.md).

[Code layout](code-layout.md) shows where each of these parts lives in the code.
