# Technical details

T-Rex is modeled from the source to the sample position, as a chain of guides with its
[`NXdisk_chopper`](../../reference/components.md) choppers in the gaps. The choppers are phased
to the central wavelength `lambda0` (see [Chopper timing](#chopper-timing)). After the sample
are a beamstop and a banana detector (see [Detectors](#detectors)). There are no polarizers or
collimators yet.
`sample_pos` is on the beam axis, 163.8 m along it from the origin.

[T-Rex](index.md) shows the whole model as McStas draws it, and [Source](source.md),
[Extraction guides](extraction.md), [Bender](bender-scan.md), [Beam monitor](beam-monitor.md)
and [Choppers](choppers.md) show what the model gives so far.
[Instrument layout](#instrument-layout) lists every component with its position. The code is in
[`mcdgs.instruments.trex`](../../reference/trex.md).

## Control parameters

| Parameter | Default | Meaning |
| --------- | ------- | ------- |
| `bender_in` | `False` | `True` uses the bender and the cold guide 4b, `False` the thermal guide 4a. |
| `b_rot` | −0.398 deg | Rotation of the bender about +y. |
| `lambda0` | 2.0 Å | Central wavelength. |
| `rrm` | 18 | Repetition rate multiplication factor. |
| `mode` | `"FF"` | `"FF"` puts the wide slits of the pulse-shaping and monochromatic choppers on the beam, for high flux; `"RR"` the narrow ones, for high resolution. |
| `t0_offset` | 1.7 ms | Time from the start of the pulse at which the neutrons the choppers are phased to leave the moderator. |
| `omega` | 0.0 deg | Goniometer rotation about +y. |

`bender_in`, `b_rot` and `omega` reach McStas as they are. `lambda0`, `rrm`, `mode` and
`t0_offset` reach it through the frequency and delay of each chopper.

## Instrument layout

`print(trex)` lists the control parameters, the McStas instrument parameters they give, and every
component in the order of the beam. For each component, `at` and `rotated` are its placement in
the `.instr` file, relative to the component named under `relative`, or to the origin where that
is `-`. `position` is where it ends up in the lab frame. Lengths are in m and angles in deg.

???+ abstract "print(trex)"

    ```text
    --8<-- "docs/instruments/trex/figures/trex_print.txt"
    ```

[`mcdgs.diagnostics.trex_00_overview`](../../reference/diagnostics.md) makes the printout, and
`pixi run docs-diagnostics` runs it.

## Source

`ESS_butterfly` on beamline W7, with its focusing rectangle at the entrance of guide 1,
1.895927 m downstream.

| Parameter | Value |
| --------- | ----- |
| Sector, beamline | W, 7 |
| Accelerator power | 2 MW |
| Moderator height | 0.03 m |
| Fraction of events from the cold moderator | 0.3 |
| Focusing rectangle | 0.094 m × 0.029992 m |
| Wavelength band | `lambda0` ± 0.905 Å, within 0.1 to 10 Å: 1.095 to 2.905 Å at the default |
| Pulses | 1 |

The source simulates only the wavelengths that can reach the sample position in one frame: those
arriving there within one source period, 1/14 Hz, plus the length of the pulse, 3.5 ms, from
t = L λ / 3956 with L = 163.8 m and λ in Å. 3956 m·Å/s is h/m_n, Planck's constant over the
neutron mass, so a neutron of 1 Å flies at 3956 m/s. The band is then

Δλ = 3956 m·Å/s × (1/14 Hz + 3.5 ms) / 163.8 m = 1.810 Å,

centered on `lambda0` and cut to 0.1 to 10 Å, the widest band of `ESS_butterfly`. The 3.5 ms is
longer than the 2.857 ms of the ESS pulse in McStas, to keep its tail. A new `lambda0`, through
`update_params`, creates the source and the beam monitor again with its band.

The source sends each ray from the moderator towards the focusing rectangle, and weights it by
the solid angle of the rectangle. The flux into guide 1 is right only if the rectangle covers the
whole entrance as seen from every point of the moderator. In the plane of the entrance it does,
when it is as large as the entrance. The rectangle is centered on the axis, while the entrance is
2.0 mm off it, from the fixed wall at x = -0.047 m to x = 0.043 m, so the rectangle is 0.094 m
wide to reach both walls.

## Extraction guides

Dimensions are from the Table of Optics V3.2 unless noted. Lengths and positions are in m, widths
and heights in mm, and z is the lab position of the entrance.

| Guide | z | Length | Width in → out | Height in → out | m left, right | m top, bottom |
| ----- | - | ------ | -------------- | --------------- | ------------- | ------------- |
| `guide_1` | 1.895927 | 0.9994 | 89.955 → 82.455 | 29.992 → 34.937 | 1.5 | 4 |
| `guide_2` | 2.896327 | 0.9992 | 82.447 → 74.948 | 34.941 → 39.160 | 1.5 | 4 |
| `guide_3` | 3.896727 | 1.4792 | 74.939 → 63.838 | 39.164 → 44.511 | 1.5 | 4 |
| `guide_4a` | 5.40065 | 0.481 | 63.838 → 60.256 | 44.594 → 46.152 | 1.5 | 3.5 |
| `guide_4b` | 5.45065 | 0.431 | 63.466 → 60.256 | 44.759 → 46.152 | 1.5 | 3.5 |

The dimensions of guide 4 are rounded here to 0.001 mm. The model uses the sheet's full values, in
`guides/extraction.py` and `guides/bispectral.py`, with each cell named in the comments.

For the top of guide 3, the sheet has 4 to 3.5 and marks it "issue". The model uses 4.

### Fixed wall

The −x wall of the NBOA (Neutron Beam Optics Assembly) is fixed at −47 mm. All of the horizontal
narrowing is taken out of the +x wall, on the cold moderator side, so the view of the thermal face
stays constant while the cold-side acceptance is trimmed.

`Guide_gravity` is symmetric about its own axis, so each guide sits on the bisector of its two
walls. The center-line is a straight line tilted by half the taper angle, not the nominal beam axis:

| Section | Taper angle | Placed on |
| ------- | ----------- | --------- |
| Guides 1 to 3 | 0.4300 deg | `arm_nboa`, tilted −0.2150 deg |
| Guide 4 | 0.4267 deg | `arm_guide_4`, tilted −0.2134 deg |

Guide 4 has its own arm because its taper differs and the width at the exit of guide 3 carries over
the gap unchanged. The arm keeps the −x wall on the same line. Guides 4a and 4b lie on one axis and
share it.

### Two paths through guide 4

Guides 4a and 4b are alternative optics in the same slot, with the same exit plane, chosen by
`bender_in`. They are not parallel beamlines.

- **Bender out** (`bender_in=0`): the thermal path, through guide 4a.
- **Bender in** (`bender_in=1`): the cold path, through the bender and guide 4b.

## Bender

`Pol_bender` at the entrance of guide 4, 0.05 m long. Its entrance is oriented toward the cold
moderator and its exit tangent is parallel to the guide axis, so nothing after it is rotated. The
nominal rotation, −0.398 deg, is the angle `atan(length / radius)`.

| Parameter | Value |
| --------- | ----- |
| Length | 0.05 m |
| Radius of curvature | 7.2 m |
| Width × height | 63.84 mm × 46 mm |
| Channels | 426 |
| Mirrors | R0 0.99, Qc 0.0219, alpha 3.02, m 4, W 0.003 |

All the mirrors have the same values for both spin states, so the bender is spin-independent.

## Curved guides

Guides 5 to 29 follow a horizontal curve of 12 km radius, from the entrance of guide 5 to the exit
of guide 29. The curve is made of straight sections, each at the angle the arc has at the section's
midpoint. `arm_start_of_curve`, at the entrance of guide 5, takes out the tilt of guide 4, since the
fixed −x wall ends there. `arm_end_of_curve`, at the exit of guide 29, turns by half the length of
guide 29, so the axis leaves the curve on its tangent.

| Guides | From, m | Width, mm | Height in → out, mm | m left, right | m top, bottom |
| ------ | ------- | --------- | ------------------- | ------------- | ------------- |
| 5 | 5.916 | 60 | 46.261 → 53.982 | 1.5 | 3.5 |
| 6 to 9 | 8.6915 | 60 | 54.015 → 76.525 | 1.5 | 3 |
| 10 | 22.7007 | 60 | 76.538 → 78.024 | 2 | 2.5 |
| 11 to 15b | 24.3017 | 60 | 78.035 → 84.965 | 1.5, 2.5 | 2.5 |
| 16 to 25 | 42.00155 | 60 | 85 | 1.5, 2.5 | 2.5 |
| 26 to 29 | 82.00155 | 60 | 85 | 1, 2.5 | 2.5 |

Positions here and in the sheet are distances along the beam axis, not lab z. The two differ by
2.8 mm at the sample position. The tables here and below summarize the sheet: the model uses its
full values, in `guides/curve.py` and `guides/straight.py`, with each cell named in the comments.

## Straight guides

The axis is straight after guide 29, at 0.4301 deg to the lab z axis. Guides 34 to 50a approximate
the ellipses of the sheet with straight tapered sections.

| Guides | From, m | Width in → out, mm | Height in → out, mm | m left, right | m top, bottom |
| ------ | ------- | ------------------ | ------------------- | ------------- | ------------- |
| 30 to 33 | 96.00255 | 60 | 85 | 1 | 2.5 |
| 34 to 37 | 111.098 | 60.129 → 69.542 | 85 | 1 | 2.5 |
| 38 to 41 | 124 | 69.542 → 66.734 | 84.998 → 78.364 | 1 | 2.5 |
| 42 | 140 | 66.734 → 63.480 | 78.364 → 74.442 | 1 | 3 |
| 43 to 45 | 144 | 63.480 → 45.385 | 74.442 → 55.335 | 1.5 | 3 |
| 47 | 155.596 | 45.385 → 38.918 | 55.335 → 49.103 | 2 | 4 |
| 48b, 49 | 158 | 38.918 → 22.433 | 49.103 → 34.987 | 2 | 4.5 |
| 50a | 162.083 | 22.433 → 14.184 | 34.987 → 29.556 | 3.5 | 4.5 |

Guide 48b is the plain guide in the polarizer slot: guides 48a and 48c are not modeled. The sheet
has no guide 46. The slit `slit_mc`, 20.5 mm × 34.5 mm, is 19 mm after guide 50a.

The arms `arm_bwc1`, `arm_bwc2`, `arm_bm1`, `arm_psc` and `arm_mc` mark the gaps for bandwidth
choppers 1 and 2, beam monitor 1, the pulse-shaping choppers and the monochromatic choppers. The
choppers `bw1`, `bw2`, `ps1`, `ps2`, `mc1` and `mc2` come right after the arm of their gap.

## Monitors

`TRex` has four monitors, `beam_monitor`, `pre_sample_monitor`, `sample_monitor` and
`beamstop_monitor`; `source_monitor` is that of [`TRexSource`](source.md). All five leave the
neutrons unchanged.
They are in [`mcdgs.instruments.trex.monitors`](../../reference/trex.md).

| Monitor | Where | Size | What it records |
| ------- | ----- | ---- | --------------- |
| `source_monitor` | The entrance of guide 1, at the place and rotation of `arm_nboa` | 0.089955 m × 0.029992 m | Time against wavelength, in bins of 0.05 ms and 0.05 Å |
| `beam_monitor` | In the middle of the bunker gap of row 29, after guide 4, 5.902 m along the beam, on the arm of guide 4 | 0.06 m × 0.0463 m | Time against wavelength, in bins of 0.05 ms and 0.02 Å |
| `pre_sample_monitor` | Between `slit_mc` and the sample, 163.5 m along the beam, on the arm of the monochromatic choppers | 0.1 m × 0.1 m | Time against wavelength, in bins of 0.05 ms and 0.02 Å |
| `sample_monitor` | At `sample_pos` | 0.1 m × 0.1 m | Time against wavelength, in bins of 0.05 ms and 0.02 Å |
| `beamstop_monitor` | On the face of the [`beamstop`](#detectors), 2.999 m after `sample_pos`, but before the sample in beam order, on the arm of the monochromatic choppers | 0.1 m × 0.1 m, that of the beamstop | Time against wavelength, in bins of 0.05 ms and 0.02 Å, of the direct beam as if there were no sample |

The time and wavelength maps (`source_monitor`, `beam_monitor`, `pre_sample_monitor`,
`sample_monitor` and `beamstop_monitor`) follow one rule. The wavelength axis is the band
of the source. The time axis runs from the fastest neutron of that band, leaving at the start of
the pulse, to the slowest, leaving 3.5 ms later. Each axis has the fewest bins of its width that
cover it, so its stop moves up to the end of the last bin. `source_monitor` has the widest band,
0.1 to 10 Å in 198 bins and 0.048 to 8.30 ms in 165. `beam_monitor` has the band around
`lambda0`: at its default, 1.095 to 2.915 Å in 91 bins and 1.63 to 7.83 ms in 124.

## Detectors

After the sample, both placed on `sample_pos`, come the beamstop, in
[`mcdgs.instruments.trex.monitors`](../../reference/trex.md) with `beamstop_monitor` on its face,
and the banana detector, in [`mcdgs.instruments.trex.detectors`](../../reference/trex.md). None
is in the Table of Optics: their sizes are placeholders, to be replaced by those of the
instrument.

In McStas the order of the components, not their place, decides what happens to a neutron.
`beamstop_monitor` comes before `sample_pos`, after `pre_sample_monitor`. Each neutron is taken
to its plane, 3 m beyond the sample, counted, and given back unchanged, still upstream of the
sample. So it counts the whole direct beam, as if there were no sample, and the sample then
scatters the same neutrons. The beamstop behind the sample stops those that pass straight
through.

| Component | Where | Size | What it does |
| --------- | ----- | ---- | ------------ |
| `beamstop` | On the beam axis, 2.999 m after `sample_pos`, 1 mm in front of the banana detector | 0.1 m × 0.1 m, as `pre_sample_monitor` | A McStas `Beamstop`: absorbs every neutron that hits it, so the direct beam does not reach the banana detector |
| `banana_detector` | A vertical cylinder of radius 3 m centered on `sample_pos` | −36° to 144° in the horizontal plane, 0° along the beam; 2.2 m high, lowered by 0.2978 m so that, seen from the sample, it spans about 25° below the beam to 15° above | Lists every neutron reaching it, with its pixel, time and wavelength. The pixels are 288 in angle by 136 in height, numbered from 0 as Mantid does |

A neutron that misses the beamstop, or is scattered away from it, goes on to the banana
detector. The beamstop is small: seen from the sample it covers ±0.95°, so the edge of the direct
beam reaches the banana detector at about ±1° to ±1.5°.

The events of `banana_detector` are read per pixel with `read_pixels` (see
[Detector pixels](../../how-to/read-the-output.md#detector-pixels)). In the tof model, the
banana detector is 3 m beyond the sample, where its neutrons arrive.

## Chopper timing

Each chopper takes its frequency and delay from the McStas instrument parameters `nu_<name>` and
`delay_<name>`, computed from `lambda0`, `rrm`, `mode` and `t0_offset`:

| Choppers | Frequency |
| -------- | --------- |
| `bw1`, `bw2` | 14 Hz, the source frequency |
| `ps1`, `ps2` | 0.75 × `rrm` × 14 Hz, counter-rotating: `ps2` turns the other way |
| `mc1`, `mc2` | `rrm` × 14 Hz, counter-rotating: they turn the same way on opposite sides of the beam |

[Choppers](choppers.md#frequencies) explains these frequencies and the counter-rotation.

A chopper's delay is the time a neutron of wavelength `lambda0` takes to reach it, plus
`t0_offset`. The flight path is the chopper's distance from the moderator along the
beam axis: the position of its gap in the Table of Optics plus its place in the gap. The path
through the bender is taken to be that through guide 4a.

The disks of the pulse-shaping and monochromatic choppers have two kinds of slit: narrow ones
(20° and 2.5°) and wide ones (35° and 4.4°). The mode chooses which is on the beam when a
neutron of wavelength `lambda0` arrives. In `"RR"` the narrow slit is centered on the beam then,
so the delay is the arrival time. In `"FF"` the disk is phased to the wide slit instead: it has
turned past its zero mark by the angle from the center of the wide slit to the beam, 55° for the
pulse-shaping choppers and 175° for the monochromatic choppers. The model keeps the disks as they
are and changes the delay, which is the same as shifting the slit edges.
[Choppers](choppers.md#delays) gives the delay of each disk at the default parameters.

[Choppers](choppers.md) draws each disk and the acceptance diagram of the cascade, and
[Simulate the choppers with tof](../../how-to/simulate-the-choppers-with-tof.md#compare-settings)
follows neutrons through it: the bandwidth choppers pass one band of wavelengths, and the
pulse-shaping and monochromatic choppers cut it into sub-pulses.

The slit edges of each disk reach McStas as an array of their own, declared in the `.instr` file
as `<name>_slit_edges`. McStas 3.10.2 sizes a list written in place by its first use, so the
2 edges of `bw1` would have been too few for the 8 of `ps1`. See the
[McStas instrument file](instr-file.md).

## Parameters file

This is the default T-Rex saved as a parameters file, `save(TRex(), PARAMS_FILE)`, the
reference of how a T-Rex file looks. `PARAMS_FILE` is `src/mcdgs/instruments/trex/params.yaml`,
and a test checks that it stays the instrument the code makes.

It holds the control parameters, then the McStas instrument parameters they give, then one
section per component in beam order: its McStas component, its parameters, where it sits (`at`,
`rotated`, `relative`) and its position in the lab frame. The instrument parameters and the
positions are for reading only, and loading the file does not read them. A file you load need
only hold what differs from this one; see [Parameters files](../../reference/io.md).
[View it on GitHub](https://github.com/bingli621/mcdgs/blob/main/src/mcdgs/instruments/trex/params.yaml),
or expand it here:

??? abstract "params.yaml"

    ```yaml
    --8<-- "src/mcdgs/instruments/trex/params.yaml"
    ```
