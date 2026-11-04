# Install and first run

This page takes you from a fresh clone to a simulated spectrum at T-Rex's sample position.

## Install

McDGS is not on PyPI. It runs from a clone of the repository, in an environment managed by
[pixi](https://pixi.sh), on Linux (x86-64) or macOS 14 or later (Apple silicon). The environment
includes McStas itself.

```bash
git clone https://github.com/bingli621/mcdgs.git
cd mcdgs
pixi install
pixi run pytest
```

!!! warning "Run everything through pixi"

    Use `pixi run python my_script.py`, or open a `pixi shell` first. Calling the environment's
    Python directly leaves out the settings that locate McStas: building an instrument stops with
    `No component named Progress_bar in McStas installation`.

## Look at T-Rex

An instrument is a Python object, and creating one does not start McStas. Printing it shows its
parameters and its components in beam order:

```python
from mcdgs.instruments.trex.model import TRex

trex = TRex()
print(trex)
```

The last table gives the placement of each component, with `position` in the lab frame:

```text
name                component      at                         rotated            relative  position
origin              Progress_bar   [0, 0, 0]                  [0, 0, 0]          -         [0, 0, 0]
source              ESS_butterfly  [0, 0, 0]                  [0, 0, 0]          -         [0, 0, 0]
arm_nboa            Arm            [-0.00202282, 0, 1.89593]  [0, -0.214998, 0]  -         [-0.00202282, 0, 1.89593]
...
sample_pos          Arm            [0, 0, 1.815]              [0, 0, 0]          arm_mc    [0.830079, 0, 163.797]
```

The whole printout is on the [T-Rex](instruments/trex/index.md#instrument-layout) page.

Each component is an attribute, and printing one shows all of its values:

```python
print(trex.guide_1)
print(trex.positions()["sample_pos"])  # [x, y, z] in m
```

## Run it

`TRex` has a wavelength monitor, `sample_monitor`, at the sample position. This script runs
T-Rex with the bender in and prints what that monitor records. Save it as `first_run.py` and run
it with `pixi run python first_run.py`.

```python
import h5py

from mcdgs.instruments.trex.model import TRex
from mcdgs.mcstas.simulation import Simulation

simulation = Simulation(TRex())
h5_file = simulation.run(ncount=1e7, bender_in=1)
print(h5_file)

with h5py.File(h5_file) as file:
    components = file["entry1/instrument/components"]
    monitor = next(c for name, c in components.items() if name.endswith("_sample_monitor"))
    spectrum = monitor["output/sample_monitor_L"]
    for wavelength, intensity in zip(spectrum["Wavelength__Angs_"][()], spectrum["data"][()]):
        print(f"{wavelength:.1f} Å  {intensity:.3g}")
```

The choppers pass one band of wavelengths around `lambda0`, 2.0 Å by default, from about 1.3 to
2.9 Å. It is highest near 2.3 Å, since the bender removes much of the short-wavelength side.
Few neutrons get through the choppers, so the script runs 1e7: with 1e6 only about a hundred
reach the monitor, and the spectrum has empty bins inside the band.

## What the run left behind

| Folder | Contents |
| ------ | -------- |
| `mcstas_input/TRex/` | The generated `.instr` and `.c` files and the compiled instrument. |
| `mcstas_output/<time>_run/` | The results of this run, in `mccode.h5`. |

Both folders are created where you ran the script, and in this repository git ignores them.

## Next

- [Core concepts](concepts.md) explains what the script above is made of.
- [Run a McStas simulation](how-to/run-a-simulation.md) covers scans and split runs.
- [Simulate the choppers with tof](how-to/simulate-the-choppers-with-tof.md) shows which
  wavelengths reach the sample, and when, in under a second.
