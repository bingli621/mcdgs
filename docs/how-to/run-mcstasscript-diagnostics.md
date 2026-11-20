# Run McStasScript's diagnostics

McStasScript has two tools to look inside a McStas run:

- its **beam diagnostics**, which record the rays at points along the instrument and draw them as
  histograms, such as the beam profile or the wavelength against time;
- its **analysis diagram**, which draws the components in beam order with the intensity and the
  number of rays just before each.

Both run McStas themselves. A [`McStasScriptDiagnostics`](../reference/mcstas.md) runs them in the
same [workspace](run-a-simulation.md#where-the-files-go) as a `Simulation`, with a run folder and a
[record](run-a-simulation.md#find-and-clean-up-runs) for every run. The examples are those of
[`mcdgs.diagnostics.trex.guides_diagnostics`](../reference/diagnostics.md), which runs `TRex` in
the white beam, 0.1 to 10 Å, with every chopper parked; [Guides](../instruments/trex/guides.md)
shows its figures and what they say about the guides of T-Rex, and the
[extraction guides](../instruments/trex/extraction.md) show both tools on the first 6 m.

```python
from mcdgs.diagnostics import cli
from mcdgs.diagnostics.trex.guides_diagnostics import WAVELENGTH_BINS, plot_transmission
from mcdgs.instruments.trex.model import Params, TRex
from mcdgs.mcstas import McStasScriptDiagnostics
```

## The beam at points along the instrument

Add the points, where the rays are recorded, and the views, how each point is drawn, then run:

<!-- fmt:off -->
```python
--8<-- "src/mcdgs/diagnostics/trex/guides_diagnostics.py:instrument"

--8<-- "src/mcdgs/diagnostics/trex/guides_diagnostics.py:diagnostics"
```
<!-- fmt:on -->

The code is that of its `main`, where `NAME` is `"trex_guides_diagnostics"`, the name of the run
folder and of the figure, `ncount` is 1e6 and `seed` is 1, and `folder` is where the figures go.
`rrm` 0 parks every chopper, `bw1` to `mc2`, with the slit of the mode on the beam: the wide
slits of `"FF"`, the default. The beam then reaches the sample position unchopped. With no
chopper to phase, it runs with `lambda0` set to `None`, the white beam, where turning choppers
would need a wavelength to be phased to.

!!! note "A chopper at 0 Hz is parked, not open"

    At 0 Hz, McStas's `NXdisk_chopper` stands still with its zero mark at `park_angle`, and a ray
    passes only through the opening in front of the beam. `rrm` 0 sets `park_<name>` to put the
    slit of the mode there. With `lambda0` at 2 Å and 1e5 neutrons, about 11000 rays reach `ps1`
    and 7200 reach `mc1`; 9 reach the sample position with the choppers turning, and about 6300
    with them parked.

- A point is just `before` or just `after` a component, by its name in the instrument:
  `list(TRex().components)` gives them all. A name the instrument does not have raises a
  `KeyError` at once.
- A view is a histogram of one ray variable, or of one against another: `x`, `y`, `z` in m;
  `vx`, `vy`, `vz` and `speed` in m/s; `dx`, `dy`, the divergence in deg; `t` in s; `l`, the
  wavelength in Å; `e`, the energy in meV. Its other keywords are those of McStasScript's
  `add_view`: `bins`; `same_scale=True` draws the view on the same ranges at every point;
  `axis1_values` and `axis2_values` draw lines at those values; and plot options, such as
  `log=True`.
- `run` takes the instrument parameters as keywords, like
  [`Simulation.run`](run-a-simulation.md#one-run), and returns the record of the run.

Then draw it, one row per point and one column per view, under the instrument's
`running_conditions()`, and save it with `save_figure`, which also closes it:

<!-- fmt:off -->
```python
--8<-- "src/mcdgs/diagnostics/trex/guides_diagnostics.py:plot"
```
<!-- fmt:on -->

`plot` draws the last run. A view of two variables is drawn in the `turbo` colormap, where
McStasScript would use `jet`; give `colormap=` to the view for another. The figure of this run is
on [Guides](../instruments/trex/guides.md#the-beam-along-the-guides).

!!! warning "Every ray is kept, so keep `ncount` low"

    A point records every ray that passes it, as McStasScript's limited lists stop the run early,
    which leaves the intensities too low. Near the source, where most rays still pass, that is
    about 60 MB per point at 1e6 neutrons, so the diagnostics run with fewer neutrons than a
    [`Simulation`](run-a-simulation.md#one-run).

!!! note "Every list of rays needs more than one ray"

    McStasScript reads back every list of rays the run writes: those of the points, and those of
    the event monitors of the instrument, such as the `banana_detector` of T-Rex. A list with
    exactly one ray fails the run as it is read. Past the turning choppers of T-Rex few rays are
    left: at 1e5 neutrons, about 9 reach `sample_pos` with the bender out, and 1 reaches the
    detector with it in. Put the points before the choppers, park them with `rrm` 0, or run
    enough neutrons for the detector.

### Histograms of the rays

`histogram` sums the rays of the last run at each point into bins of one ray variable, in n/s per
bin, or counts them with `weighted=False`. It gives a dict by point, such as `"after source"` or
`"before bw1"`, in the order of the instrument. `guides_diagnostics` divides the wavelength
spectra of two points for the transmission of the guides between them:

<!-- fmt:off -->
```python
--8<-- "src/mcdgs/diagnostics/trex/guides_diagnostics.py:transmission"
```
<!-- fmt:on -->

The figure is on [Guides](../instruments/trex/guides.md#by-wavelength).

`rays` gives the rays themselves: for each point, a dict of the ray variables asked for, an
array each with a value per ray, empty at a point no ray reached. Besides the variables of the
views, `"p"` is the weight in n/s, and `"dx"` and `"dy"` the horizontal and vertical divergence
in deg. [`brilliance_transfer`](../reference/diagnostics.md) reads them before the beam monitor
and before `sample_pos`, and keeps the rays in a window of position and divergence:

<!-- fmt:off -->
```python
--8<-- "src/mcdgs/diagnostics/trex/brilliance_transfer.py:trex"
```
<!-- fmt:on -->

The figures are on [Brilliance transfer](../instruments/trex/brilliance-transfer.md).

### The same diagnostics on more runs

The points and the views stay, so one setup runs again with other instrument parameters. Every run
builds the instrument anew and writes its own folder:

```python
diag.run("bender_out", ncount=1e6, seed=1, bender_in=0)
bender_out = diag.plot()
diag.run("bender_in", ncount=1e6, seed=1, bender_in=1)
bender_in = diag.plot()
```

For a change that is not an instrument parameter, such as another component, make a new
instrument and a new `McStasScriptDiagnostics`, as the
[extraction guides](../instruments/trex/extraction.md) do on each path.

## The intensity before each component

`run_analysis` runs McStas once more, with a monitor before every component, and draws
McStasScript's diagram of the components with the intensity, in n/s, and the number of rays before
each, on a log scale. Each step of the line is what one component loses:

<!-- fmt:off -->
```python
--8<-- "src/mcdgs/diagnostics/trex/guides_diagnostics.py:analysis"
```
<!-- fmt:on -->

It needs no points or views. The diagram of this run is on
[Guides](../instruments/trex/guides.md#the-intensity-before-each-component), which reads the
transmission of each section of guide from it.

The arrows on the left go from each component to the one it is placed relative to. A component
with a dashed box takes part only under a condition, such as guide 4a, the bender and guide 4b
of T-Rex, which depend on `bender_in`; the line is not to be read there.

McStasScript has no range for the intensity
axis, so `intensity_limits`, in n/s, sets it, such as to see the losses along the guides:

```python
run, fig = diag.run_analysis(ncount=1e7, seed=1, intensity_limits=(1e10, 2e12))
```

The ray count keeps its ratio to the intensity, so both share their ticks.

The monitors before the components are 100 m wide, so they count every ray, inside an opening or
not. The folder of the run also holds the monitors of the instrument.

## Where the runs go

Each run writes a folder in `runs/` of the workspace, like every run of a `Simulation`:
`<time>_diagnose` for `run` and `<time>_analysis` for `run_analysis`, unless given a name, such as
`trex_guides_diagnostics` above. A `folder_name` is replaced on every run. They are found and
cleaned up by their type:

```python
from mcdgs.mcstas import RunType

workspace = diag.workspace
workspace.latest_run(RunType.DIAGNOSE)  # the newest run with the diagnostics
workspace.clean_runs(keep=2, run_type=RunType.ANALYSIS)
```

Both write McStas's own format, which McStasScript reads, not McStas NeXus: their record names a
`data_file` that does not exist. Read the monitors of a run with McStasScript instead:

```python
import mcstasscript as ms

data = ms.load_data(run.folder)
```

The rays of the last `run` are in `diag.diagnostics`, McStasScript's own `Diagnostics`.
