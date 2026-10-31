# Add an instrument

An instrument is a package in `src/mcdgs/instruments/`, with a control-parameters dataclass and a
subclass of `Instrument`. This small example is complete and runs.

## Declare the control parameters

`param` gives each field a default, a unit and a comment, which `print(instrument)` shows.

```python
from dataclasses import dataclass

from mcdgs.core.params import param


@dataclass(frozen=True)
class Params:
    """Control parameters of the instrument."""

    width: float = param(0.1, "m", "width of the detector")
```

## Create the components

Set `params_class` and `name`, then create the components in `__init__`, in beam order, after
calling `super().__init__(params)`. McStas places a component relative to earlier ones, so a
`relative` must name a component above it.

`__init__` must accept `params`. `instr_params` says which parameters a run may set.

```python
from mcdgs.core.base import Instrument
from mcdgs.core.components import Arm, ESS_butterfly, Monitor_nD, Progress_bar
from mcdgs.core.params import InstrParam


class Small(Instrument):
    """A source, a sample position and a detector."""

    params_class = Params
    name = "small"  # the McStas instrument name

    def __init__(self, params=None):
        super().__init__(params)
        self.origin = Progress_bar()
        self.source = ESS_butterfly(Lmin=2, Lmax=4, dist=2, focus_xw=0.05, focus_yh=0.05)
        self.sample_pos = Arm(at=[0, 0, 2], rotated=[0, 0.5, 0])
        self.detector = Monitor_nD(
            xwidth="width",  # an instrument parameter, as in McStas
            yheight=0.1,
            options="x bins=4, y bins=3, list all, pixel id, t",
            at=[0, 0, 1],
            relative="sample_pos",
        )

    @property
    def instr_params(self):
        return [InstrParam("width", self.params.width, "m", "width of the detector")]
```

A component cannot take a name the instrument already uses, such as `name` or `positions`: that
raises a `ValueError`.

## Check it

Print it to see the positions, then run a few neutrons:

```python
from mcdgs.mcstas.simulation import Simulation

small = Small()
print(small)
h5_file = Simulation(small).run(ncount=1e4, width=0.05)
```

## A component that does not exist yet

Add a dataclass to the matching module of `src/mcdgs/core/components/`, named after the McStas
component, and list it in `__init__.py` there. Copy `Arm` or `Monitor_nD`:

- use `@dataclass(kw_only=True, repr=False)` and inherit from `Component`;
- name the fields after the McStas parameters, with McStas's defaults;
- leave out parameters you do not need: McStas keeps its defaults.

## Finish

- Add tests in `tests/instruments/`, as in `tests/instruments/trex/test_trex.py`.
- Add a page in `docs/instruments/` and a page in `docs/reference/`, and list them in the `nav` of
  `zensical.toml`.
