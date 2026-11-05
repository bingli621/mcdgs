# Source

`TRexSource`, in
[`mcdgs.diagnostics.trex.source`](https://github.com/bingli621/mcdgs/blob/main/src/mcdgs/diagnostics/trex/source.py),
compares McStas and [tof](../../how-to/simulate-the-choppers-with-tof.md) at the entrance of
guide 1, the first optics of T-Rex. Nothing stands between the source and that point, so the two
should see the same pulse there: what differs is only how each models the ESS moderator.

## The setup

It is the smallest T-Rex: two components and no McStas instrument parameters.

| Component | What it is |
| --------- | ---------- |
| `source` | The T-Rex source, `ESS_butterfly` on beam port W7, from 0.1 to 10 Å, one pulse at 2 MW, the same component as in `TRex`. |
| `source_monitor` | A `Monitor_nD` at the entrance of guide 1, 1.895927 m from the source on its axis, unrotated, with the opening of guide 1: 0.089955 m wide and 0.029992 m high. |

The monitor records one 2D map of time against wavelength:

| Axis | Range | Bins | Bin width |
| ---- | ----- | ---- | --------- |
| Time | 0 to 15 ms | 300 | 0.05 ms |
| Wavelength | 0 to 10 Å | 200 | 0.05 Å |

15 ms holds the whole pulse and the flight to guide 1 up to 10 Å. The monitor has
`restore_neutron=1`, so it does not change the beam.

`TRexSource` has the control parameters `lambda0` and `rrm`, with the defaults of `TRex`. Nothing
uses them yet: they set no McStas parameter and no component.

It is defined like any instrument of mcdgs: a parameters class, then an `Instrument` subclass that
creates its components in beam order, with their positions.

<!-- fmt:off -->
```python
--8<-- "src/mcdgs/diagnostics/trex/source.py:trex_source"
```
<!-- fmt:on -->

## Running both models

The same `TRexSource` object runs in McStas and in tof. McStas writes the instrument file, runs it
and saves the monitor in its run folder, which McStasScript reads back as n/s per bin:

<!-- fmt:off -->
```python
--8<-- "src/mcdgs/diagnostics/trex/source.py:run_mcstas"
```
<!-- fmt:on -->

tof builds its model from the same components: the `ESS_butterfly` becomes tof's pulse of beam
port W7 and the `Monitor_nD` a detector named `source_monitor`. It gives each neutron's time and
wavelength there, which `read_tof_monitor` counts on the bins of the McStas monitor:

<!-- fmt:off -->
```python
--8<-- "src/mcdgs/diagnostics/trex/source.py:run_tof"
```
<!-- fmt:on -->

Both runs take the same `seed`. See [Run a McStas simulation](../../how-to/run-a-simulation.md)
and [Simulate the choppers with tof](../../how-to/simulate-the-choppers-with-tof.md) for each.

### Where the tof source sits

The two models place their source differently:

| | Source at | Monitor at | Flight to the monitor |
| --- | --- | --- | --- |
| McStas | 0 m, the moderator | 1.895927 m | 1.895927 m |
| tof | 0.05 m | 1.895927 m | 1.845927 m |

McStas emits each neutron at the moderator, with its time of emission. tof's pulse of beam port
W7, `ess-trex`, was not sampled there. It was sampled 5 cm away, so its times are those of the
neutrons when they pass 0.05 m, and they already include the flight over the first 5 cm. The
builder therefore keeps the tof source at 0.05 m, and the monitor at its distance from the
moderator, 1.895927 m. Both models then have the neutrons reach the monitor at the same time.

With the source at 0 m, tof would add the flight over those 5 cm a second time. That is
25 µs at 2 Å and 126 µs at 10 Å, half a bin to two and a half bins of 0.05 ms in the figure
below.

Every other distance in tof is the straight line from the moderator to the component. For
T-Rex that is within 1 mm of the path along the curved guide. Up to guide 1 the path is straight,
so here they are the same.

## Counts and intensity

Each simulated neutron stands for part of the real flux, its weight. The intensity of a bin is its
counts times the weight, in n/s. The counts depend only on how many neutrons are simulated; the
intensity is the physics, and should not.

- **McStas** sets the weights itself, so that the rays add up to the flux of the moderator. With
  ten times the `ncount`, a bin has ten times the rays, each with a tenth of the weight: the
  intensity stays the same, and only its noise goes down.
- **tof** draws its neutrons from a distribution with no flux: each has weight 1, so its
  intensity would be its counts, and would grow with the number of neutrons drawn. tof therefore
  takes its weight from McStas:

    ```text
    scale = McStas total (n/s) / tof total (counts)
    ```

    is the n/s that each tof neutron stands for. Scaled, tof has the level of McStas, and the shape
    in time and wavelength of its own model.

The comparison is thus of shapes, not of absolute flux: by construction both have the same total.

## Error bars

A bin of N neutrons has the Poisson error √N, a relative error of 1/√N. In n/s it is √N times the
weight. With a hundred times the neutrons, the error in counts is ten times larger, but the
relative error is a tenth, and so is the error in n/s, as the weight is a hundredth. The error
bars shrink with more neutrons in both models, while the intensity stays the same.

McStas gives the error of each bin itself. For tof it is √N of the counts, scaled by the same
factor as the counts. The bins are independent, so when the maps are summed over one axis their
errors add in quadrature.

## The figure

`main` runs both models with the same seed, reads each monitor on the same bins, prints the mean
time and wavelength of each and the McStas total, scales tof to McStas, and saves
`trex_source.png` next to the `mccode.h5` of the McStas run. The figure has two rows:

- **Top:** the time and wavelength map of McStas and of tof, side by side on one color scale, in
  n/s per 0.05 ms and 0.05 Å.
- **Bottom:** the wavelength spectrum, summed over time, and the time profile, summed over
  wavelength, of both models on the same axes, with their error bars. The legend gives the
  McStas `ncount`, the number of tof neutrons and the n/s that each of them stands for.

![McStas and tof at the entrance of guide 1: time and wavelength maps, the wavelength spectrum and the time profile](figures/trex_source.png)

This is `main` at its defaults: 1e7 McStas neutrons, 1,000,000 tof neutrons and seed 1. The total
at the entrance of guide 1 is 8.24e11 n/s, so each tof neutron stands for 8.24e5 n/s. Both models
have the same means, 2.84 ms and 2.40 Å.

The two models agree within their error bars almost everywhere. The spectrum has the thermal
peak near 1.1 Å, the cold one near 2.6 Å, and a third peak at the shortest wavelengths, below
0.3 Å. Only there do they differ: in the first bins McStas is about 10 % higher than tof. In time,
both rise over the 2.86 ms of the pulse and peak near 2.8 ms; the long wavelengths arrive later,
up to about 8 ms at 10 Å.

The figure needs McStas, so it is not made with the site: `pixi run docs-mcstas` makes it again,
with the other figures of T-Rex that need McStas.

To make it yourself, run:

```bash
pixi run python -m mcdgs.diagnostics.trex.source
```

This uses 1e6 McStas and 10,000 tof neutrons, enough for a quick look. `main` defaults to 1e7
and 1,000,000; call it with other numbers for a finer figure:

```bash
pixi run python -c "from mcdgs.diagnostics.trex.source import main; main(ncount=1e8, neutrons=10_000_000)"
```
