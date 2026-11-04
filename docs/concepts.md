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

## Three kinds of parameter

| Kind | Where it lives | Changed by | Example |
| ---- | -------------- | ---------- | ------- |
| Design parameter | The fields of a component | Editing the instrument | `guide_1.w1` |
| Control parameter | `instrument.params` | `instrument.update_params(...)` | `bender_in=True` |
| McStas instrument parameter | `instrument.instr_params`, computed from the control parameters | A keyword of a run | `simulation.run(bender_in=1)` |

The instrument parameters are the ones McStas knows about. Their values in `instr_params` are the
defaults of a run, and a keyword of a run overrides them for that run only, without changing the
instrument.

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
   and a detector at every monitor and at each component named in `detectors`. It leaves out the
   optics, and measures each distance in a straight line from the source.
2. `model.run()` follows the neutrons in memory and writes no files.
3. `result["sample_pos"]` holds the neutrons at that detector, with their wavelength and time of
   arrival. See [Simulate the choppers with tof](how-to/simulate-the-choppers-with-tof.md).

## Layers

Dependencies point inwards, so an instrument never depends on how it is simulated:

| Package | Role |
| ------- | ---- |
| `core` | `Instrument`, the components, parameters and geometry. It imports nothing else from `mcdgs`. |
| `instruments`, `samples` | The instruments and ready-made crystals, built on `core`. |
| `mcstas` | An adapter that works with any instrument: builds it, runs it, reads the output. |
| `tof` | An adapter that works with any instrument: builds its chopper cascade as a tof model. |
| `io` | An adapter that saves any instrument to a YAML parameters file and loads it back. |
| `reduction` | Reserved for data reduction. It holds no code yet. |
| `diagnostics` | Setups that study one part of an instrument, such as T-Rex up to guide 4, and what reads their McStas and tof results. One package per instrument, built on `instruments`; nothing else imports it. |
