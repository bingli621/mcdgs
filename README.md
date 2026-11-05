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
pixi run docs-mcstas            # remake the files of the site that need McStas, then commit them
```

The site makes what needs only mcdgs and tof, such as the table of the chopper modes, from the
current code for every build. What needs McStas, such as the McStas drawings of T-Rex and the
figures of its diagnostics, is committed: `docs-mcstas` runs every diagnostics module that needs
McStas and writes their files to `docs/instruments/trex/figures/`.

## Layout

```
src/mcdgs/
  core/          entities shared by any instrument: Instrument, components, params, geometry
  instruments/   the instruments, one package each
  samples/       ready-made samples
  io/            instruments to and from dicts and YAML files
  mcstas/        build the McStasScript instrument, run simulations, read NeXus output
  tof/           build the chopper cascade of an instrument as a tof model
  reduction/     data reduction
  diagnostics/   setups that study one part of an instrument, one package per instrument; each
                 module with a main makes the files of one page of the documentation
scripts/         build tooling of the documentation site
notebooks/       Jupyter notebooks of one-off simulations that diagnose the instrument models
tests/
docs/            the documentation site, built with Zensical
mcstas_input/    generated McStas files, per instrument (git-ignored)
mcstas_output/   simulation results, a folder per run (git-ignored)
```

Dependencies point inwards: `core` imports nothing from the other layers; `instruments` and
`samples` build on it; `io`, `mcstas`, `tof` and `reduction` are adapters that work with any
instrument. `diagnostics` builds on `instruments` and reads McStas and tof results, and nothing in
`mcdgs` imports it.

## License

[BSD-3-Clause](LICENSE)
