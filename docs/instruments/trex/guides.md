# Guide transmission

[`mcdgs.diagnostics.trex.guides_diagnostics`](https://github.com/bingli621/mcdgs/blob/main/src/mcdgs/diagnostics/trex/guides_diagnostics.py)
follows the beam of T-Rex through all its guides, from the source to the sample position. It
runs `TRex` in the white beam, 0.1 to 10 Å (`lambda0` `None`), with every chopper parked on its
wide slit (`rrm` 0, mode `"FF"`), so the beam is not chopped: what it loses on the way is what
the guides take out, and the little the parked disks block. The bender is out, so the beam goes
through the thermal guide 4a.

It uses McStasScript's two tools: its beam diagnostics, the beam at five points along the
instrument, and its diagram of the intensity before every component. From the rays of the
diagnostics it also draws the transmission of the guides by wavelength. How much of the
brilliance of the moderator the guides keep, in windows of position and divergence at the
sample position, and on both paths through guide 4, is on
[Brilliance transfer](brilliance-transfer.md).
[Run McStasScript's diagnostics](../../how-to/run-mcstasscript-diagnostics.md) has the code and
explains both tools; this page shows what they give for T-Rex. The guides themselves are listed
in [Technical details](details.md#curved-guides).

## The beam along the guides

The points are after the source, before the first bandwidth chopper `bw1` (32.0 m), after guide
29 at the end of the curve (96.0 m), after guide 50a at the end of the straight guide (163.1 m)
and before `sample_pos` (163.8 m). Each row is one point; the columns are the beam profile,
wavelength against time, and the spectra in time and in wavelength, in n/s per bin. At 1e6
neutrons and seed 1:

![The beam after the source, before the first bandwidth chopper, at the end of the curved and of the straight guide, and at the sample position](figures/guides/trex_guides_diagnostics.png)

After the source, the beam leaves the two moderators side by side, brighter at x above 0,
within the pulse, about 3 ms long. Its spectrum has two peaks, at about 1.1 and 2.5 Å, and a
long tail to 10 Å. Further along, a neutron arrives later the longer its wavelength, so
wavelength against time is a band. The band keeps the length of the pulse while the flight time
grows with the distance, so it looks thinner at every point. Before the first bandwidth chopper,
the beam fills the guide, ±0.03 m by ±0.04 m, and has 14 times less intensity than after the
source. From there to the sample position, the spectrum keeps its two peaks, and drops by a
factor of about 2.3; 114224 of the 1e6 rays reach the sample position. There, the beam is a spot
about 0.04 m across.

## The intensity before each component

McStasScript's diagram lists the components in beam order, with the intensity just before each,
in n/s, and the number of rays, on a log scale. Each step of the line is what one component
loses. At 1e7 neutrons and seed 1:

??? abstract "The diagram of the intensity before each component"

    ![The components of T-Rex and the intensity before each](figures/guides/trex_intensity.png)

The intensity falls from 8.8e11 n/s at the entrance of guide 1 to 2.5e10 n/s at the sample
position. Most of it is lost in the first guides, up to about 20 m; the parked choppers take
little. Guide 4a, the bender and guide 4b have dashed boxes: they take part only under a
condition on `bender_in`, and the line is not to be read there. At `arm_bwc1`, `arm_bwc2` and
`arm_end_of_curve` the line dips to about half and comes back at the next component. An arm
takes nothing from the beam, so the dip is how the diagram's monitor reads at those arms, not a
loss. With the choppers turning, the pulse-shaping and monochromatic choppers take the largest
steps, down to about 1e7 n/s at the detector.

## Transmission of the guides

### Along the beam

The transmission of a section is the intensity at its exit over that at its entrance, read from
the diagram's run, at 1e7 neutrons. The sections follow the tables of
[Technical details](details.md#curved-guides):

| Guides | From → to, m | In, n/s | Out, n/s | Transmission |
| ------ | ------------ | ------- | -------- | ------------ |
| 1 to 4a, extraction | 1.9 → 5.9 | 8.84e11 | 1.48e11 | 0.167 |
| 5 to 9, start of the curve | 5.9 → 22.7 | 1.48e11 | 5.97e10 | 0.404 |
| 10 to 25, curve, with `bw1` and `bw2` | 22.7 → 82.0 | 5.97e10 | 5.16e10 | 0.863 |
| 26 to 29, end of the curve | 82.0 → 96.0 | 5.16e10 | 4.39e10 | 0.851 |
| 30 to 33, straight, with `ps1` and `ps2` | 96.0 → 111.1 | 4.39e10 | 3.73e10 | 0.849 |
| 34 to 45, elliptic | 111.1 → 155.6 | 3.73e10 | 3.54e10 | 0.950 |
| 47 to 50a, focusing, with `mc1` and `mc2` | 155.6 → 163.1 | 3.54e10 | 2.48e10 | 0.701 |
| All | 1.9 → 163.8 | 8.84e11 | 2.48e10 | 0.028 |

From the beam monitor after guide 4, the guides pass 0.168 of the beam to the sample position.
The losses are not spread evenly:

- **The first 20 m take most.** The extraction guides pass a sixth of the beam, and guides 5 to
  9 less than half of the rest, most of it in guides 5 and 6. Most likely, the beam leaves the
  extraction with more divergence than guides 5 to 9, with m 1.5 on their sides, carry.
- **The middle of the curve loses little**, under 1 % per 4 m guide. Guides 26 to 29 lose more,
  where the m of one side drops from 1.5 to 1, and so do guides 30 to 33, the straight guide
  after the curve, with m 1 on both sides: they trim the divergence the curve leaves.
- **Guides 34 to 45 pass 95 %.** Most of what reaches them already fits within m 1.
- **The focusing end takes 30 %**, most in guide 49, 17.5 %, and guide 50a, 12 %: the guides
  narrow from 45 × 55 mm to 14 × 30 mm, which raises the divergence of the beam.
- **The parked choppers take about 3 % together**: 0.1 % each for `bw1` and `bw2`, 1.1 % for the
  pulse-shaping pair and 1.7 % for the monochromatic pair. Parked on their wide slits, they
  block little of the beam.

### By wavelength

The rays of the beam diagnostics give the transmission by wavelength, as the ratio of the
intensities at two points, in bins of 0.1 Å, at 1e6 neutrons. On the left are the sections
between the points, on the right the whole way to the sample position, on a log scale. A bin
that fewer than 20 rays reach has no value, which leaves out most wavelengths below 0.4 Å:

![The transmission by wavelength of the sections of guide, and from the source and from bw1 to the sample position](figures/guides/trex_guides_transmission.png)

The end of the straight guide, after guide 50a, is left out: the beam reaches the sample
position from there unchanged. In wider bands, from the same rays:

| From → to | 0.1–0.5 Å | 0.5–1 Å | 1–1.5 Å | 1.5–2 Å | 2–3 Å | 3–4 Å | 4–6 Å | 6–8 Å | 8–10 Å |
| --------- | --------- | ------- | ------- | ------- | ----- | ----- | ----- | ----- | ------ |
| Source → `bw1` | 0.009 | 0.027 | 0.050 | 0.061 | 0.058 | 0.078 | 0.122 | 0.207 | 0.300 |
| `bw1` → end of the curve | 0.151 | 0.808 | 0.886 | 0.891 | 0.825 | 0.727 | 0.674 | 0.646 | 0.646 |
| End of the curve → sample | 0.324 | 0.487 | 0.637 | 0.627 | 0.601 | 0.554 | 0.523 | 0.492 | 0.457 |
| `bw1` → sample | 0.049 | 0.393 | 0.564 | 0.559 | 0.496 | 0.403 | 0.353 | 0.318 | 0.295 |

Here the source is the beam just after it, all the source sends toward guide 1. The mean
wavelength is 2.60 Å after the source, 3.58 Å before `bw1`, 3.36 Å at the end of the curve and
3.20 Å at the sample position.

- **Up to `bw1`, the long wavelengths pass best**, from under 1 % below 0.5 Å to 30 % at 8 to
  10 Å. A supermirror of a given m reflects up to an angle in proportion to the wavelength, so
  the first guides carry the long wavelengths at much larger angles than the short ones.
- **After `bw1`, 1 to 2 Å pass best**, 56 % to the sample position, and the transmission falls
  toward long wavelengths, to 30 % at 8 to 10 Å. The long wavelengths the first guides took in
  at large angles reflect more often, and lose more at each reflection. Gravity adds to it: a
  neutron of 8 Å, at 494 m/s, takes 0.27 s from `bw1` to the sample and would fall 0.35 m, four
  times the height of the guide, so it bounces off the bottom wall all the way.
- **Below 0.5 Å, the curve takes almost everything.** Over its 90 m the 12 km curve leaves no
  line of sight through the 60 mm guide, so a neutron must reflect off its walls to pass. The
  angle that needs, about √(2 × 60 mm / 12 km) = 0.18°, is the critical angle of m 2.5 at
  about 0.7 Å: shorter wavelengths fall off, to 15 % below 0.5 Å.

Overall, from the source to the sample position, about 3 % pass at 1 to 4 Å, both peaks
alike, and 9 % at 8 to 10 Å: the first guides favor the long wavelengths more than the rest
favor the short ones, and the mean wavelength rises from 2.60 Å to 3.20 Å.

## Making the figures

The figures need McStas, so they are not made with the site: `pixi run docs-diagnostics` makes
them again, with the other files of the T-Rex pages. To make them yourself, run:

```bash
pixi run python -m mcdgs.diagnostics.trex.guides_diagnostics
```

It saves the three figures to `mcdgs_workspace/results/trex_guides_diagnostics/`, or to the
folder given with `--out`. The transmissions of the tables are read from the run folders it
writes, `trex_guides_diagnostics` and `trex_guides_diagnostics_analysis` in the workspace.
