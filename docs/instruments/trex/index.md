# T-Rex

T-Rex is modelled from the source to the sample position, as a chain of guides. There are no
choppers, polarizers, collimators or detector yet, though the chopper gaps are in place as arms
and the [`NXdisk_chopper`](../../reference/components.md) component is ready for them.
`sample_pos` is on the beam axis, 163.8 m along it from the origin.

[Extraction results](extraction.md) shows what the model gives so far. The code is in
[`mcdgs.instruments.trex`](../../reference/trex.md).

## Control parameters

| Parameter | Default | Meaning |
| --------- | ------- | ------- |
| `bender_in` | `False` | `True` uses the bender and the cold guide 4b, `False` the thermal guide 4a. |
| `b_rot` | −0.398 deg | Rotation of the bender about +y. |
| `lambda0` | 2.0 Å | Central wavelength. |
| `rrm` | 18 | Repetition rate multiplication factor. |
| `omega` | 0.0 deg | Goniometer rotation about +y. |

Only `bender_in` and `b_rot` reach McStas. Nothing in the model uses `lambda0`, `rrm` or `omega`
yet.

## Source

`ESS_butterfly` on beamline W7, with a focusing rectangle 2 m downstream.

| Parameter | Value |
| --------- | ----- |
| Sector, beamline | W, 7 |
| Accelerator power | 2 MW |
| Moderator height | 0.03 m |
| Fraction of events from the cold moderator | 0.3 |
| Focusing rectangle | 0.095 m × 0.035 m |
| Wavelength band | 0.1 to 10 Å |
| Pulses | 1 |

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
`guides/extraction.py`, with each cell named in the comments.

For the top of guide 3, the sheet has 4 to 3.5 and marks it "issue". The model uses 4.

### Fixed wall

The −x wall of the NBOA (Neutron Beam Optics Assembly) is fixed at −47 mm. All of the horizontal
narrowing is taken out of the +x wall, on the cold moderator side, so the view of the thermal face
stays constant while the cold-side acceptance is trimmed.

`Guide_gravity` is symmetric about its own axis, so each guide sits on the bisector of its two
walls. The centre-line is a straight line tilted by half the taper angle, not the nominal beam axis:

| Section | Taper angle | Placed on |
| ------- | ----------- | --------- |
| Guides 1 to 3 | 0.4300 deg | `arm_NBOA`, tilted −0.2150 deg |
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

| Parameter | Value | Source |
| --------- | ----- | ------ |
| Length | 0.05 m | Table of Optics |
| Radius of curvature | 7.2 m | `trex_primary.txt` |
| Width × height | 63.84 mm × 46 mm | `trex_guide_anyshape.ipynb` |
| Channels | 426 | `trex_primary.txt` |
| Mirrors | R0 0.99, Qc 0.0219, alpha 3.02, m 4, W 0.003 | `trex_primary.txt` |

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
2.8 mm at the sample position. The tables here and below summarise the sheet: the model uses its
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

Guide 48b is the plain guide in the polarizer slot: guides 48a and 48c are not modelled. The sheet
has no guide 46. The slit `slit_mc`, 20.5 mm × 34.5 mm, is 19 mm after guide 50a. Its size is from
`trex_primary.txt`.

The arms `arm_bwc1`, `arm_bwc2`, `arm_bm1`, `arm_psc` and `arm_mc` mark the gaps for bandwidth
choppers 1 and 2, beam monitor 1, the P-Chopper and the M-Chopper.
