# McDGS

McDGS is a Monte Carlo (MC) simulation of the Direct Geometry Spectrometer (DGS) at the European
Spallation Source (ESS), from [McStas](https://www.mcstas.org) to S(Q,E). An instrument is
described once, as a Python class holding its components and parameters. From that description,
McDGS builds the McStas instrument and runs the simulation through
[McStasScript](https://github.com/PaNOSC-ViNYL/McStasScript), then loads the detector data from
the McStas NeXus files with [McStasToX](https://github.com/mads-bertelsen-agentic/McStasToX). The
reduction will normalize the data by the incident monitor and the pixel solid angle and convert
it to momentum and energy transfer, with algorithms based on [Scipp](https://scipp.github.io),
[ScippNeutron](https://scipp.github.io/scippneutron) and [tof](https://scipp.github.io/tof).

```python
from mcdgs.instruments.trex.model import TRex

trex = TRex()
print(trex)  # every parameter and component, with its position in the lab frame
```

## Status

| Part | State |
| ---- | ----- |
| [T-Rex](instruments/trex/index.md) | Modelled from the source to the sample position, guides only. The choppers, the polarizers, the collimators and the detector are not modelled yet. |
| Simulation | Single runs, parameter scans and split runs at the sample position. |
| Reading the output | Monitors and detector pixels from McStas NeXus files. |
| Instrument files | Not written yet. |
| Reduction to S(Q,E) | Not written yet. |

## Where to go next

- [Install and first run](getting-started.md): set up the environment and run T-Rex's extraction
  guides.
- [Core concepts](concepts.md): how instruments, components and parameters fit together.
- [T-Rex](instruments/trex/index.md) and its
  [extraction results](instruments/trex/extraction.md): the model and what it shows so far.
- [How-to guides](how-to/run-a-simulation.md): run a scan, read the output, add an instrument.
- [Reference](reference/core.md): the Python API, generated from the source.
