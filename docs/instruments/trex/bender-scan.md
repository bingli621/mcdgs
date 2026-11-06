# Bender scan

[`mcdgs.diagnostics.trex_03_bender_scan`](https://github.com/bingli621/mcdgs/blob/main/src/mcdgs/diagnostics/trex_03_bender_scan.py)
scans the rotation of the [bender](details.md#bender), `b_rot`, and shows how it passes each band
of wavelengths. It runs the [`Extraction`](extraction.md) and reads its beam monitor, after
guide 4.

## The scan

`run_bender_scan` runs McStas at 19 values of `b_rot`, from −0.8 to 0.1 deg in steps of 0.05 deg,
with the bender in, and once with it out, through the thermal guide 4a. `read_spectra` sums the
map of the beam monitor over time, and then over bands of 0.5 Å centered on 0.5, 1.0, … 9.5 Å:
19 bands over the 0.1 to 10 Å of the monitor. Each wavelength bin, 0.02 Å wide, goes to the band
whose center is nearest to its own, 25 bins to a band. A bin is centered on each edge of two
bands, and goes to the band above, so each band is that of its center moved up by half a bin,
0.01 Å. The bins are independent, so their errors add in quadrature.
`transmission` divides each band of the scan by the same band with the bender out, and the
relative errors of the two add in quadrature.

![Intensity at the beam monitor against the rotation of the bender, per wavelength band, over the intensity with the bender out](figures/bender_scan.png)

Each curve is one band, its color the center of the band on the colorbar. The dashed line at 1 is
the bender-out intensity, and the dotted line the nominal rotation, −0.398 deg.

## Transmission at the nominal rotation

At −0.4 deg, the point of the scan nearest the nominal rotation:

| Band center | Transmission |
| ----------- | ------------ |
| 8.5 to 9.5 Å | 0.81 to 0.86 |
| 3.0 to 8.0 Å | 0.93 to 1.02 |
| 2.5 Å | 0.88 |
| 2.0 Å | 0.68 |
| 1.5 Å | 0.29 |
| 0.5 and 1.0 Å | 0.13 and 0.12 |

The errors are 0.007 to 0.04, the larger at the longer wavelengths, where fewer neutrons are.
Above 1 means the bender passes a band slightly better than guide 4a; the bands from 3 to 8 Å
pass within a few percent of guide 4a, within about two errors. They are flat from about −0.6 to −0.2 deg, so the nominal rotation
sits in the middle of that range. Above 8 Å the transmission is lower, 0.81 to 0.86 here and
0.75 to 0.93 over the scan. The bands of 1.0 and 1.5 Å keep rising up to the end of the scan at
0.1 deg, so their maximum lies outside it. The shortest band, around 0.5 Å, reaches at most 0.2.

## Making the figure

The figure needs McStas, so it is not made with the site: `pixi run docs-diagnostics` makes it
again, with the other files of the T-Rex pages. To make it yourself, run:

```bash
pixi run python -m mcdgs.diagnostics.trex_03_bender_scan
```

This is `main` at its defaults: 1e6 McStas neutrons per angle, and seed 1. The scan goes to
`mcstas_output/bender_scan/`, the run without the bender to `mcstas_output/bender_out/`, each
replacing the one of an earlier run, and the figure, `bender_scan.png`, next to the scan.
`read_angles`, `read_spectra` and `transmission` read the two files again, without running McStas.
