# McDGS

McDGS is a Monte Carlo (MC) simulation of the Direct Geometry Spectrometer (DGS) at the European
Spallation Source (ESS), from [McStas](https://www.mcstas.org) to S(Q,E). An instrument is
described once, as a Python class holding its components and parameters. That one description
can be simulated in two ways:

- **McStas**, the whole instrument: McDGS builds the McStas instrument and runs it through
  [McStasScript](https://github.com/PaNOSC-ViNYL/McStasScript), then loads the detector data from
  the McStas NeXus files with [McStasToX](https://github.com/mads-bertelsen-agentic/McStasToX).
- **[tof](https://scipp.github.io/tof)**, the chopper cascade alone: the source, the choppers and
  the flight path, without the optics, in under a second.

The reduction will normalize the data by the incident monitor and the pixel solid angle and convert
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
| [T-Rex](instruments/trex/details.md) | Modeled from the source to the sample position: the guides and the choppers. The polarizers, the collimators and the detector are not modeled yet. |
| Simulation with McStas | Single runs, parameter scans and split runs at the sample position. |
| [Simulation with tof](how-to/simulate-the-choppers-with-tof.md) | The source, the choppers and any monitor or named position, with the choppers timed as in McStas. |
| Reading the output | Monitors and detector pixels from McStas NeXus files; the neutrons at each tof detector. |
| [Parameters files](reference/io.md) | Save an instrument to a YAML file and load it back. A file need only hold what differs from the instrument the code makes. |
| Reduction to S(Q,E) | Not written yet. |

## Where to go next

- [Install and first run](getting-started.md): set up the environment and run T-Rex to the sample
  position.
- [Core concepts](concepts.md): how instruments, components and parameters fit together.
- [T-Rex](instruments/trex/details.md) and its
  [extraction results](instruments/trex/extraction.md): the model and what it shows so far.
- [How-to guides](how-to/run-a-simulation.md): run a McStas scan, simulate the choppers with tof,
  read the output, add an instrument.
- [Reference](reference/core.md): the Python API, generated from the source.
