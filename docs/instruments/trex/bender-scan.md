# Bender scan

[`mcdgs.diagnostics.trex.bender_scan`](https://github.com/bingli621/mcdgs/blob/main/src/mcdgs/diagnostics/trex/bender_scan.py)
scans the rotation of the [bender](details.md#bender), `b_rot`, and shows how it passes each band
of wavelengths. It runs the [`Extraction`](extraction.md) and reads its beam monitor, after
guide 4.

## The scan

`run_bender_scan` runs McStas at 19 values of `b_rot`, from −0.8 to 0.1 deg in steps of 0.05 deg,
with the bender in, and once with it out, through the thermal guide 4a. `read_spectra` sums the
map of the beam monitor over time, and then over bands of 8 of its wavelength bins: 12 bands of
0.442 Å, from 0 to 5.3 Å. The bins are independent, so their errors add in quadrature.
`transmission` divides each band of the scan by the same band with the bender out, and the
relative errors of the two add in quadrature.

![Intensity at the beam monitor against the rotation of the bender, per wavelength band, over the intensity with the bender out](figures/bender_scan.png)

Each curve is one band, labeled by its center. The dashed line at 1 is the bender-out intensity,
and the dotted line the nominal rotation, −0.398 deg. The colors repeat after ten bands: the
lower blue and orange curves are 0.22 and 0.66 Å, the upper ones 4.64 and 5.08 Å.

## Transmission at the nominal rotation

At −0.4 deg, the point of the scan nearest the nominal rotation:

| Band center | Transmission |
| ----------- | ------------ |
| 2.87 Å and above | 0.94 to 1.06 |
| 2.43 Å | 0.88 |
| 1.99 Å | 0.67 |
| 1.55 Å | 0.33 |
| 0.22 to 1.10 Å | 0.11 to 0.14 |

The errors are 0.01 to 0.05, the larger at the longer wavelengths, where fewer neutrons are.
Above 1 means the bender passes a band slightly better than guide 4a; within the errors the bands
of 2.9 Å and above pass completely. They are flat from about −0.6 to −0.2 deg, so the nominal
rotation sits in the middle of that range. The bands from 0.66 to 1.55 Å keep rising up to the end
of the scan at 0.1 deg, so their maximum lies outside it. The shortest band, around 0.22 Å, stays
near 0.1 at every angle.

## Making the figure

The figure needs McStas, so it is not made with the site: `pixi run docs-diagnostics` makes it
again, with the other files of the T-Rex pages. To make it yourself, run:

```bash
pixi run python -m mcdgs.diagnostics.trex.bender_scan
```

This is `main` at its defaults: 1e6 McStas neutrons per angle, and seed 1. The scan goes to
`mcstas_output/bender_scan/`, the run without the bender to `mcstas_output/bender_out/`, each
replacing the one of an earlier run, and the figure, `bender_scan.png`, next to the scan.
`read_angles`, `read_spectra` and `transmission` read the two files again, without running McStas.
