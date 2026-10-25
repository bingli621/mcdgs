# mcdgs

[![CI](https://github.com/bingli621/mcdgs/actions/workflows/ci.yml/badge.svg)](https://github.com/bingli621/mcdgs/actions/workflows/ci.yml)

## Use

```bash
pixi run pytest --cov=mcdgs     # tests and coverage
pixi run pre-commit install     # once: ruff format on every commit
```

## Layout

```
src/mcdgs/
  core/          entities shared by any instrument: Instrument, components, params, geometry
  instruments/   the instruments, one package each
  samples/       ready-made samples
  io/            instruments to and from dicts and YAML files
  mcstas/        build the McStasScript instrument, run simulations, read NeXus output
  reduction/     data reduction
scripts/         build, run and plot
tests/
mcstas_input/    generated McStas files, per instrument (git-ignored)
mcstas_output/   simulation results, a folder per run (git-ignored)
output/          figures
```

Dependencies point inwards: `core` imports nothing from the other layers; `instruments` and
`samples` build on it; `io`, `mcstas` and `reduction` are adapters that work with any
instrument.
