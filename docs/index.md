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

The reduction normalizes the data by the incident monitor and the pixel solid angle and converts
it to momentum and energy transfer, with the kernels of [Scipp](https://scipp.github.io) and
[ScippNeutron](https://scipp.github.io/scippneutron).

```python
from mcdgs.instruments.trex.model import TRex

trex = TRex()
print(trex)  # every parameter and component, with its position in the lab frame
```

## Status

| Part | State |
| ---- | ----- |
| [T-Rex](instruments/trex/details.md) | Modeled from the source to the detector: the guides, the choppers, the monitors, a beamstop and a banana detector. The polarizers and the collimators are not modeled yet. |
| Simulation with McStas | Single runs, parameter scans and split runs at the sample position. |
| [Simulation with tof](how-to/simulate-the-choppers-with-tof.md) | The source, the choppers and every monitor, with the choppers timed as in McStas. |
| [Chopper acceptance](instruments/trex/choppers.md#acceptance-diagram) | The acceptance diagram of the choppers, exact and without neutrons, with scippneutron. |
| Reading the output | Monitors and detector pixels from McStas NeXus files; the neutrons at each tof detector. |
| [Parameters files](reference/io.md) | Save an instrument to a YAML file and load it back. A file need only hold what differs from the instrument the code makes. |
| [Reduction](reference/reduction.md) | The incident energy of each sub-pulse from two monitors, and the energy and momentum transfer of each neutron, per incident neutron and per steradian: the [elastic plane](instruments/trex/elastic-plane.md) of a crystal. S(Q,E) is not written yet. |

## Where to go next

- [Install and first run](getting-started.md): set up the environment and run T-Rex to the sample
  position.
- [Core concepts](concepts.md): how instruments, components and parameters fit together.
- [T-Rex](instruments/trex/details.md) and its results so far: the
  [source](instruments/trex/source.md), the [extraction guides](instruments/trex/extraction.md),
  the [bender](instruments/trex/bender-scan.md), the
  [beam monitor](instruments/trex/beam-monitor.md), the [choppers](instruments/trex/choppers.md),
  the [neutrons on the sample](instruments/trex/neutrons-on-sample.md), the
  [masked source](instruments/trex/masked-source.md), the
  [incident energy](instruments/trex/incident-energy.md), and the
  [elastic plane](instruments/trex/elastic-plane.md) and
  [rotation scan](instruments/trex/elastic-scan.md) of a crystal.
- [How-to guides](how-to/run-a-simulation.md): run a McStas scan, simulate the choppers with tof,
  read the output.
- [Reference](reference/core.md): the Python API, generated from the source.
