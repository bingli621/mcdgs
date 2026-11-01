# Install and first run

This page takes you from a fresh clone to a simulated spectrum at the exit of T-Rex's guide 4.

## Install

mcdgs is not on PyPI. It runs from a clone of the repository, in an environment managed by
[pixi](https://pixi.sh), on Linux (x86-64) or macOS (Apple silicon). The environment includes
McStas itself.

```bash
git clone https://github.com/bingli621/mcdgs.git
cd mcdgs
pixi install
pixi run pytest
```

!!! warning "Run everything through pixi"

    Use `pixi run python my_script.py`, or open a `pixi shell` first. McStas compiles each
    instrument, and calling the environment's Python directly leaves out the settings the compiler
    needs: the run stops with `'napi.h' file not found`.

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
name         component      at                         rotated            relative     position
origin       Progress_bar   [0, 0, 0]                  [0, 0, 0]          -            [0, 0, 0]
source       ESS_butterfly  [0, 0, 0]                  [0, 0, 0]          -            [0, 0, 0]
arm_NBOA     Arm            [-0.00202282, 0, 1.89593]  [0, -0.214998, 0]  -            [-0.00202282, 0, 1.89593]
guide_1      Guide_gravity  [0, 0, 0]                  [0, 0, 0]          arm_NBOA     [-0.00202282, 0, 1.89593]
...
```

Each component is an attribute, and printing one shows all of its values:

```python
print(trex.guide_1)
print(trex.positions()["guide_4a"])  # [x, y, z] in m
```

## Run it

`TRex` has no monitor yet, so a run would record nothing. This script defines the same beamline up
to guide 4 with a wavelength monitor at the exit, then runs it with the bender in. Save it as
`first_run.py` and run it with `pixi run python first_run.py`.

```python
import h5py

from mcdgs.core.base import Instrument
from mcdgs.core.components import Monitor_nD, Progress_bar
from mcdgs.core.params import InstrParam
from mcdgs.instruments.trex.guides.extraction import guides_extraction
from mcdgs.instruments.trex.model import Params, source
from mcdgs.mcstas.simulation import Simulation


class Extraction(Instrument):
    """T-Rex up to guide 4, with a wavelength monitor at the exit."""

    params_class = Params
    name = "my_extraction"

    def __init__(self, params=None):
        super().__init__(params)
        self.origin = Progress_bar(at=[0, 0, 0])
        self.source = source()
        for name, component in guides_extraction().items():
            setattr(self, name, component)
        self.spectrum = Monitor_nD(
            xwidth=0.1,
            yheight=0.1,
            restore_neutron=1,
            options="lambda limits=[0.5 5.3] bins=12",
            filename="spectrum",
            at=[0, 0, self.guide_4b.l + 1e-3],
            relative="guide_4b",
        )

    @property
    def instr_params(self):
        return [
            InstrParam("bender_in", int(self.params.bender_in), "", "bender in = 1, out = 0"),
            InstrParam("b_rot", self.params.b_rot, "deg", "rotation of the cold neutron bender"),
        ]


simulation = Simulation(Extraction())
h5_file = simulation.run(ncount=1e6, bender_in=1)
print(h5_file)

with h5py.File(h5_file) as file:
    spectrum = file["entry1/instrument/components/0010_spectrum/output/spectrum_L"]
    for wavelength, intensity in zip(spectrum["Wavelength__Angs_"][()], spectrum["data"][()]):
        print(f"{wavelength:.1f} Å  {intensity:.3g}")
```

McStas's compiler warnings about static vectors come from the bender and are harmless. The
spectrum peaks near 2.7 Å, the cold peak: the bender has removed most of the short wavelengths.

## What the run left behind

| Folder | Contents |
| ------ | -------- |
| `mcstas_input/my_extraction/` | The generated `.instr` and `.c` files and the compiled instrument. |
| `mcstas_output/<time>_run/` | The results of this run, in `mccode.h5`. |

Both folders are created where you ran the script, and in this repository git ignores them.

## Next

- [Core concepts](concepts.md) explains what the script above is made of.
- [Run a simulation](how-to/run-a-simulation.md) covers scans and split runs.
