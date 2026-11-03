# McDGS

[![CI](https://github.com/bingli621/mcdgs/actions/workflows/ci.yml/badge.svg)](https://github.com/bingli621/mcdgs/actions/workflows/ci.yml)

Monte Carlo (MC) simulation of the Direct Geometry Spectrometer (DGS) at the European Spallation
Source. From McStas to S(Q,E).

**Documentation: <https://bingli621.github.io/mcdgs/>**, including the T-Rex model and its results.

## Use

```bash
pixi run pytest --cov=mcdgs     # tests and coverage
pixi run pre-commit install     # once: ruff format on every commit
pixi run docs                   # preview the documentation site
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
scripts/         build tooling of the documentation site
notebooks/       Jupyter notebooks of one-off simulations that diagnose the instrument models
tests/
docs/            the documentation site, built with Zensical
mcstas_input/    generated McStas files, per instrument (git-ignored)
mcstas_output/   simulation results, a folder per run (git-ignored)
```

Dependencies point inwards: `core` imports nothing from the other layers; `instruments` and
`samples` build on it; `io`, `mcstas` and `reduction` are adapters that work with any
instrument.

## License

[BSD-3-Clause](LICENSE)
