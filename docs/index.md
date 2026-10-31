# mcdgs

Monte Carlo simulation of the Direct Geometry Spectrometer (DGS) at the European Spallation
Source (ESS). From McStas to S(Q,E).

An instrument is described once, as a Python class holding its components and parameters.
mcdgs builds the McStas instrument from that description, runs it and reads the output.

```python
from mcdgs.instruments.trex.trex import TRex

trex = TRex()
print(trex)  # every parameter and component, with its position in the lab frame
```

## Status

| Part | State |
| ---- | ----- |
| [T-Rex](instruments/trex/index.md) | Modelled from the source to guide 4, including the bender. The rest of the instrument is not modelled yet. |
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
