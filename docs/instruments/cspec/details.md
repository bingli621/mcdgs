# Technical details

CSPEC is modeled from the source to the sample position, as a chain of guides with its
[`NXdisk_chopper`](../../reference/components.md) choppers in the gaps. The choppers are phased
to the central wavelength `lambda0` (see [Chopper timing](#chopper-timing)). There are no
monitors, detectors, beamstop or guide exchanger yet. `sample_pos` is on the trajectory of the
last guide, 160.0 m from the moderator along it.

The dimensions are from CSPEC_nl-specs_1.E.xlsx unless noted. Numbers that are not in it are
marked FIXME in the code, and listed under
[Not in CSPEC_nl-specs_1.E.xlsx](#not-in-the-sheet).

[CSPEC](index.md) shows the whole model as McStas draws it, and [Source](source.md) shows the
beam at the entrance of the guide system. [Instrument layout](#instrument-layout) lists every
component with its position. The code is in
[`mcdgs.instruments.cspec`](../../reference/cspec.md).

## Control parameters

| Parameter | Default | Meaning |
| --------- | ------- | ------- |
| `t0_offset` | 0.0 s | Time from the start of the pulse at which the neutrons the choppers are phased to leave the moderator. |
| `lambda0` | 3.0 Å | Central wavelength. |
| `nu_mc` | 14 Hz | Frequency of the monochromatic choppers. |
| `n` | 1 | The RRM chopper turns at `nu_mc / n`. |
| `omega` | 0.0 deg | Goniometer rotation about +y. |

`omega` reaches McStas as it is. `lambda0`, `nu_mc`, `n` and `t0_offset` reach it through the
frequency and delay of each chopper.

## Instrument layout

`print(cspec)` lists the control parameters, the McStas instrument parameters they give, and
every component in the order of the beam. For each component, `at` and `rotated` are its
placement in the `.instr` file, relative to the component named under `relative`, or to the
origin where that is `-`. `position` is where it ends up in the lab frame. Lengths are in m and
angles in deg.

??? abstract "print(cspec)"

    ```text
    --8<-- "docs/instruments/cspec/figures/overview/cspec_print.txt"
    ```

[`mcdgs.diagnostics.cspec.overview`](../../reference/diagnostics.md) makes the printout, and
`pixi run docs-diagnostics` runs it.

## Source

`ESS_butterfly` on beamline W3, at 138°, with its focusing rectangle in the plane of the entrance
of the guide system, 1.90424 m downstream.

| Parameter | Value |
| --------- | ----- |
| Sector, beamline | W, 3 |
| Accelerator power | 2 MW |
| Moderator height | 0.03 m |
| Fraction of events from the cold moderator | 0.9 (FIXME) |
| Focusing rectangle | 0.0776 m × 0.043 m |
| Wavelength band | `lambda0` ± 0.927 Å, within 0.1 to 10 Å: 2.074 to 3.926 Å at the default |
| Pulses | 1 |

The source simulates only the wavelengths that can reach the sample position in one frame: those
arriving there within one source period, 1/14 Hz, plus the length of the pulse, 3.5 ms, from
t = L λ / 3956 with L = 160.004 m and λ in Å. 3956 m·Å/s is h/m_n, Planck's constant over the
neutron mass. The band is then

Δλ = 3956 m·Å/s × (1/14 Hz + 3.5 ms) / 160.004 m = 1.853 Å,

centered on `lambda0` and cut to 0.1 to 10 Å, the widest band of `ESS_butterfly`. A new
`lambda0`, through `update_params`, gives the source its band.

The entrance of the guide system, the opening of W03-01-011, is 45 mm wide and 43 mm high, and
16.3 mm off the axis of the source. The focusing rectangle is centered on the axis, so it is
2 × (16.3 + 22.5) = 77.6 mm wide, to reach the far wall of the entrance on both sides, and as
high as the entrance. Seen from every point of the moderator it then covers the whole opening,
and the flux into the guide is right.

CSPEC_nl-specs_1.E.xlsx does not give the fraction of events from the cold moderator. At 0.9 the
spectrum of McStas is that of tof's pulse of beam port W3, with its thermal peak: see
[Source](source.md#over-the-widest-band).

## The ISCS frame

CSPEC_nl-specs_1.E.xlsx gives each element of the guide in ISCS, the frame of the insert of beam
port W3, in mm: x along the insert, y across it, z up. Each element has its start in ISCS and
two angles, those of the trajectory at its middle: χ, down and back up in the vertical
S-bender, and φ, sideways in the horizontal curve. The element runs straight from its start at
those angles.

`arm_iscs` is that frame in the McStas frame of the source. It is at the origin of ISCS,
(0.032710, 0, 0.024373) m from the focal point of the sector, and turned by −0.5° about y: the
insert is at 137.5°, the beamline at 138°. In `arm_iscs`, McStas x is ISCS y, McStas y is ISCS z
and McStas z is ISCS x, and an element is `ROTATED (−χ, −φ, 0)`. Every guide, and every arm of a
chopper gap, is placed relative to `arm_iscs` from the start and the angles of its row. The
placements agree with the McStas columns of CSPEC_nl-specs_1.E.xlsx to 22 µm, which round the ISCS
origin to 0.1 mm.

The gaps of 1 to 40 mm between the elements are left empty.

## Guides

Every guide is a `Guide_gravity`, tapered linearly from its entrance to its exit. Positions are
lab z of the entrance, in m; widths and heights in mm. The tables summarize CSPEC_nl-specs_1.E.xlsx:
the model uses their full values, in `guides/bender.py` and `guides/straight.py`, with each cell
named in the comments.

### S-bender

From the entrance of the guide system to 54.4 m the trajectory turns down and back up in a
vertical S-bender of 1600 m radius, which takes the beam 422 mm down, and sideways in a
horizontal curve of 4000 m radius, by 0.745°.

| Guides | From, m | Width in → out | Height in → out | m left, right | m top, bottom |
| ------ | ------- | -------------- | --------------- | ------------- | ------------- |
| W03-01-01 to 01-02 | 1.904 | 45 → 70 | 43 → 55 | 3, 3.5 | 4, 3.5 |
| W03-01-03 to 06-011 | 3.904 | 70 | 55 | 3, 3.5 | 3.5, 3 |
| W03-06-012 to 06-13 | 28.402 | 70 | 55 | 2.5, 3 | 3, 3.5 |

W03-01-01 is one element in two rows, and CSPEC_nl-specs_1.E.xlsx gives only its entrance and its
exit. The model splits it into `guide_01_011` and `guide_01_012`, with the middle where a
linear taper puts it (FIXME).

### Straight guides

After the S-bender the trajectory is straight, 422 mm below the beamline and turned 0.745° from
the insert. The guides approximate the ellipses of CSPEC_nl-specs_1.E.xlsx with straight tapered
elements.

| Guides | From, m | Width in → out | Height in → out | m left, right | m top, bottom |
| ------ | ------- | -------------- | --------------- | ------------- | ------------- |
| W03-06-14 to 06-331 | 54.402 | 70 → 114.66 | 55 → 111.28 | 2 | 2 |
| W03-06-332 to 07-01 | 93.593 | 114.66 | 111.28 → 73.89 | 2.5 | 2.5 |
| W03-08-01 to 08-11 | 105.760 | 114.66 → 106.33 | 73.89 → 73.6 | 2 | 2.5 to 3 |
| W03-08-12 to 08-26 | 128.736 | 106.33 → 27.54 | 73.6 → 49.86 | 2 to 3.5 | 3.5 |
| W03-09-01, 09-02 | 157.975 | 27.44 → 22.58 | 49.84 → 49.0 | 3.5 | 3.5 |

The guide exchanger after W03-09-02 is not modeled: CSPEC_nl-specs_1.E.xlsx leaves it to be defined.
The beam flies free over the 1.2 m from the exit of W03-09-02 to `sample_pos`.

## Choppers

The arms `arm_bwc1`, `arm_bwc2`, `arm_bwc3`, `arm_psc` and `arm_mc` mark the start of the gaps
of the bandwidth choppers 1 to 3, the pulse-shaping choppers and the monochromatic choppers.
Each chopper comes right after the arm of its gap, its distance from the start of the gap along
it. Its flight path is the trajectory from the moderator to the start of the gap, from
CSPEC_nl-specs_1.E.xlsx, plus that distance.

| Chopper | Gap | Into the gap, mm | Flight path, m | Slits | Axle | Beam window |
| ------- | --- | ---------------- | -------------- | ----- | ---- | ----------- |
| `bw1` | `arm_bwc1` | 15 | 14.93342 | 1 of 40.7° | Above | 70 mm, W03-05-01 |
| `bw2` | `arm_bwc2` | 15 | 20.45545 | 1 of 41.7° | Above | 70 mm, W03-05-05 |
| `bw3` | `arm_bwc3` | 15 | 104.61213 | 1 of 193° | Above | 114.66 mm, W03-06-39 |
| `ps1` | `arm_psc` | 12.71 | 105.66872 | 3 of 24.23°, 120° apart | Above | 114.66 mm, W03-07-01 |
| `ps2` | `arm_psc` | 22.71 | 105.67872 | 3 of 24.23°, 120° apart | Below | 114.66 mm, W03-07-01 |
| `rrm` | `arm_mc` | 24.11 | 158.47993 | 1 of 4.45° | Below | 22.43 mm, W03-09-01 |
| `mc1` | `arm_mc` | 48.22 | 158.50404 | 1 of 4.45° | Above | 22.43 mm, W03-09-01 |
| `mc2` | `arm_mc` | 58.22 | 158.51404 | 1 of 4.45° | Below | 22.43 mm, W03-09-01 |

CSPEC_nl-specs_1.E.xlsx gives where the first disk of each chopper turns: `bw1`, `bw2`, `bw3`, `ps1`
and `mc1`. The rest is not in it:

- The second disk of each counter-rotating pair, `ps2` and `mc2`, is 10 mm behind the first
  (FIXME).
- The RRM chopper turns halfway from the start of the gap of the monochromatic choppers to
  `mc1` (FIXME). At the start itself it would be 27 nm before the exit of W03-09-01: McStas would
  take it to be behind the neutrons leaving the guide, and absorb them all.
- The slits and the axles, above or below the beam, are those of the CSPEC model of drcspec
  (FIXME). The disks are 0.35 m in radius, with slits 0.1 m high (FIXME).

Each chopper has a beam window as wide as the exit of the guide in front of it, on the source
side.

## Chopper timing

Each chopper takes its frequency and delay from the McStas instrument parameters `nu_<name>` and
`delay_<name>`, computed from `lambda0`, `nu_mc`, `n` and `t0_offset`, as in drcspec:

| Choppers | Frequency |
| -------- | --------- |
| `bw1`, `bw2`, `bw3` | 14 Hz, the source frequency |
| `ps1`, `ps2` | `nu_mc` / 2, counter-rotating: `ps2` turns the other way (FIXME) |
| `rrm` | `nu_mc` / `n` |
| `mc1`, `mc2` | `nu_mc`, counter-rotating: `mc2` turns the other way |

Each disk has the center of its first slit on the beam when a neutron of wavelength `lambda0`,
leaving the moderator at `t0_offset`, reaches it. A disk whose axle is above the beam has the
beam at 180° from its zero mark, one whose axle is below at 0°. The delay is the time the
neutron takes to reach the disk, plus `t0_offset`, less the time the disk takes to turn from its
zero mark to the first slit:

delay = `t0_offset` + L / v(`lambda0`) − phase / (360° × ν),

with L the flight path, ν the signed frequency and phase the angle from the center of the first
slit to the beam, 180° for a disk with its axle above and 0° below.

The slit edges of each disk reach McStas as an array of their own, declared in the `.instr`
file as `<name>_slit_edges`. See the [McStas instrument file](instr-file.md).

## Not in CSPEC_nl-specs_1.E.xlsx { #not-in-the-sheet }

| What | Value in the model |
| ---- | ------------------ |
| Fraction of events from the cold moderator | 0.9 |
| Middle of W03-01-01 | Where a linear taper puts it |
| Second disk of the pulse-shaping and monochromatic pairs | 10 mm behind the first |
| Place of the RRM chopper | Halfway from the start of its gap to `mc1` |
| Slits and axles of the disks | Those of drcspec |
| Radius of the disks, height of the slits | 0.35 m, 0.1 m |
| Frequency of the pulse-shaping choppers | Half of `nu_mc` |
| Guide exchanger | Not modeled |

## Parameters file

This is the default CSPEC saved as a parameters file, `save(CSPEC(), PARAMS_FILE)`, the
reference of how a CSPEC file looks. `PARAMS_FILE` is `src/mcdgs/instruments/cspec/params.yaml`,
and a test checks that it stays the instrument the code makes.

It holds the control parameters, then the McStas instrument parameters they give, then one
section per component in beam order: its McStas component, its parameters, where it sits (`at`,
`rotated`, `relative`) and its position in the lab frame. The instrument parameters and the
positions are for reading only, and loading the file does not read them. A file you load need
only hold what differs from this one; see [Parameters files](../../reference/io.md).
[View it on GitHub](https://github.com/bingli621/mcdgs/blob/main/src/mcdgs/instruments/cspec/params.yaml),
or expand it here:

??? abstract "params.yaml"

    ```yaml
    --8<-- "src/mcdgs/instruments/cspec/params.yaml"
    ```
