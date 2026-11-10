# McDGS

[![CI](https://github.com/bingli621/mcdgs/actions/workflows/ci.yml/badge.svg)](https://github.com/bingli621/mcdgs/actions/workflows/ci.yml)

Monte Carlo (MC) simulation of the Direct Geometry Spectrometer (DGS) at the European Spallation
Source. From McStas to S(Q,E).

**Documentation: <https://bingli621.github.io/mcdgs/>**, including the T-Rex model and its results.

## Use

```bash
pixi run pytest --cov=mcdgs     # tests and coverage
pixi run ruff check             # lint; CI also runs ruff format --check
pixi run mypy src/mcdgs         # types
pixi run pre-commit install     # once: ruff format on every commit
pixi run docs                   # preview the documentation site
pixi run docs-diagnostics       # remake the files of the site, then commit them
```

The files the site includes, such as the table of the chopper modes, the McStas drawings of
T-Rex and the figures of its diagnostics, are committed: the site does not make them.
`docs-diagnostics` runs every diagnostics module and writes their files to
`docs/instruments/trex/figures/`.

## Layout

```
src/mcdgs/
  core/          entities shared by any instrument: Instrument, components, params, geometry
  instruments/   the instruments, one package each: T-Rex so far, its source, guides,
                 choppers and monitors
  samples/       ready-made samples
  io/            instruments to and from dicts and YAML files
  mcstas/        build the McStasScript instrument, run simulations, read NeXus output
  tof/           the chopper cascade of an instrument: a tof model, scippneutron's disk
                 choppers and the acceptance diagram
  reduction/     data reduction with scippneutron: the incident energy of the sub-pulses,
                 the energy and momentum transfer of each neutron, the pixel solid angle
  diagnostics/   setups that study one part of an instrument, one module each, numbered per
                 instrument (trex_01_source, ...); each module with a main makes the files of
                 one page of the documentation
scripts/         build tooling of the documentation site, and the pixi activation that
                 registers the plugins of graphviz
tests/
docs/            the documentation site, built with Zensical
mcstas_input/    generated McStas files, per instrument (git-ignored)
mcstas_output/   simulation results, a folder per run (git-ignored)
```

Dependencies point inwards: `core` imports nothing from the other layers; `instruments` and
`samples` build on it; `io`, `mcstas`, `tof` and `reduction` are adapters that work with any
instrument. `io` imports `instruments` only so that a file can name any of them, and `reduction`
imports nothing else from `mcdgs`. `diagnostics` builds on all of them, and nothing in `mcdgs`
imports it.

## License

[BSD-3-Clause](LICENSE)
