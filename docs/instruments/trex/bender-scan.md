# Bender

[`mcdgs.diagnostics.trex.bender_scan`](https://github.com/bingli621/mcdgs/blob/main/src/mcdgs/diagnostics/trex/bender_scan.py)
shows what the [bender](details.md#bender) does to the beam. Its instrument, `TRexBenderScan`, is
T-Rex up to the beam monitor after guide 4, as the [`Extraction`](extraction.md), with the bender
in, through the bender and the cold guide 4b. Only its beam monitor differs: it covers 0.25 to
9.75 Å, the bands the scan reads, instead of 0.1 to 10 Å. It scans the rotation of the bender, `b_rot`, and
reads the beam monitor, to show how it passes each band of wavelengths. The path through the
bender at its nominal rotation, below, is run by [Extraction guides](extraction.md).

<!-- fmt:off -->
```python
--8<-- "src/mcdgs/diagnostics/trex/bender_scan.py:instrument"
```
<!-- fmt:on -->

## Bender in at the nominal rotation

With the bender in at its nominal rotation, the beam goes through the bender and the cold guide
4b. `main` of [Extraction guides](extraction.md) runs this path after the one without the bender,
with the same runs of McStas and tof: what is said there of the axes, the monitors and the scale
of tof holds here too. These figures have 1e6 McStas neutrons in the run with the
diagnostics, 1e7 in the run of the diagram and the monitors, 1e6 tof neutrons and seed 1.

<!-- fmt:off -->
```python
--8<-- "src/mcdgs/diagnostics/trex/extraction.py:paths"
```
<!-- fmt:on -->

### Beam along the guides

The rows are the beam after the source, after guide 3 and at the exit of guide 4b; the columns
the beam profile, wavelength against x, and the spectrum.

![Beam after the source, guide 3 and guide 4 with the bender in](figures/extraction/trex_extraction_bender_in_diagnose.png)

The top two rows are those without the bender. At the exit of guide 4b the thermal peak is mostly
gone and the cold peak stays.

### Intensity along the guides

![The components and the intensity before each, with the bender in](figures/extraction/trex_extraction_bender_in_intensity.png)

| Lost in | Share | Intensity after it |
| ------- | ----- | ------------------ |
| guide 1 | 48 % | 4.56e11 n/s |
| guide 2 | 36 % | 2.90e11 n/s |
| guide 3 | 41 % | 1.71e11 n/s |
| bender | 36 % | 1.11e11 n/s |
| guide 4b | 6 % | 1.04e11 n/s |

Guides 1 to 3 are as without the bender. Guide 4a is skipped, so its row is flat. The bender takes
36 % and guide 4b 6 %, against the 14 % of guide 4a, and 1.04e11 n/s reaches the beam monitor,
12 % of the beam before guide 1 and of that through the entrance of guide 1, against 17 % without
the bender.

### McStas and tof at the beam monitor

1.03e11 n/s reach the monitor, against 1.47e11 n/s without the bender. tof is scaled as without
the bender, with a `transmission` of 0.2.

![McStas and tof at the beam monitor with the bender in](figures/extraction/trex_extraction_bender_in.png)

At the cold peak, 2.6 Å, McStas and tof agree, as without the bender, and above about 4 Å McStas
is above tof. McStas has almost nothing below about 1.5 Å, and its time profile peaks about
1.5 ms later than tof's, since the beam is made of longer, slower wavelengths.

## The scan

`run_bender_scan` runs McStas at 19 values of `b_rot`, from −0.8 to 0.1 deg in steps of 0.05 deg,
with the bender in, and once with it out, through the thermal guide 4a. The scan is one
`Simulation.scan` of the `Extraction`, with `bender_in` and `b_rot` as McStas instrument
parameters, so McStas compiles it once:

<!-- fmt:off -->
```python
--8<-- "src/mcdgs/diagnostics/trex/bender_scan.py:run_scan"
```
<!-- fmt:on -->

`plot_bender_scan` reads the map of the beam monitor of every scan point into scipp, along
`b_rot`, and that of the run without the bender, both with `read_mcstas_maps` (see [Read the output](../../how-to/read-the-output.md#into-scipp)). It sums each map over time,
and then over bands of 0.5 Å centered on 0.5, 1.0, … 9.5 Å: 19 bands, from 0.25 to 9.75 Å. The
beam monitor covers just these, in bins of 0.02 Å, so the edge of every band is the edge of a
bin, and each band sums 25 whole bins. The bins are independent, so their errors add in
quadrature, as scipp adds their variances:

<!-- fmt:off -->
```python
--8<-- "src/mcdgs/diagnostics/trex/bender_scan.py:read"
```
<!-- fmt:on -->

Each band of the scan is then divided by the same band with the bender out, and the
relative errors of the two add in quadrature. scipp does not divide every point by the one
reference with its errors, since that would correlate the points, so each point gets a copy of it
of its own:

<!-- fmt:off -->
```python
--8<-- "src/mcdgs/diagnostics/trex/bender_scan.py:transmission"
```
<!-- fmt:on -->

![Intensity at the beam monitor against the rotation of the bender, per wavelength band, over the intensity with the bender out](figures/bender-scan/bender_scan.png)

Each curve is one band, named in the legend by its center and half width. The first ten bands
are solid lines and the next ones dash-dotted, in the same colors again. The dashed line at 1 is
the bender-out intensity, and the dotted line the nominal rotation, −0.398 deg.

### Transmission at the nominal rotation

At −0.4 deg, the point of the scan nearest the nominal rotation:

| Band center | Transmission |
| ----------- | ------------ |
| 8.5 to 9.5 Å | 0.83 to 0.87 |
| 3.0 to 8.0 Å | 0.88 to 1.01 |
| 2.5 Å | 0.89 |
| 2.0 Å | 0.65 |
| 1.5 Å | 0.29 |
| 0.5 and 1.0 Å | 0.16 and 0.13 |

The errors are 0.007 to 0.04, the larger at the longer wavelengths, where fewer neutrons are.
Above 1 means the bender passes a band slightly better than guide 4a; the bands from 3 to 8 Å
pass 0.95 of guide 4a on average, each within a few errors of it. Their average is flat, 0.94 to
0.96, from −0.6 to −0.2 deg, so the nominal rotation sits in the middle of that range. Above 8 Å
the transmission is lower, 0.83 to 0.87 here and 0.76 to 0.96 over the scan. The bands of 1.0 and 1.5 Å keep rising up to the end of the scan at
0.1 deg, so their maximum lies outside it. The shortest band, around 0.5 Å, reaches at most 0.19.

## Making the figures

The figures need McStas, so they are not made with the site: `pixi run docs-diagnostics` makes
them again, with the other files of the T-Rex pages. To make them yourself, run:

```bash
pixi run python -m mcdgs.diagnostics.trex.bender_scan
```

This is `main` at its defaults: 1e6 McStas neutrons per angle, and seed 1. The scan goes to
`mcdgs_workspace/runs/bender_scan/`, the run without the bender to
`mcdgs_workspace/runs/bender_out/`, each replacing the one of an earlier run, and the figure,
`bender_scan.png`, to `mcdgs_workspace/results/trex_bender_scan/`.
`plot_bender_scan` reads the two runs again, without running McStas: give it
`define_instrument()` and the records of the runs, `Run.load("mcdgs_workspace/runs/bender_scan")`
and the same of `bender_out`.
The figures of the path with the bender in are made by
[Extraction guides](extraction.md#making-the-figures).
