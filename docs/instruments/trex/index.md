# T-Rex

T-Rex is modelled from the source up to guide 4. The rest of the instrument is not modelled yet:
there are no choppers, no guides after guide 4 and no detector. `sample_pos` marks the current
estimate of the sample position, 163.8 m along z.

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
`guides_extraction.py`, with each cell named in the comments.

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
| Width × height | 63.84 mm × 46 mm | `trex_guide_anyshape.ipynb` of mctrex |
| Channels | 426 | `trex_primary.txt` |
| Mirrors | R0 0.99, Qc 0.0219, alpha 3.02, m 4, W 0.003 | `trex_primary.txt` |

All the mirrors have the same values for both spin states, so the bender is spin-independent.
