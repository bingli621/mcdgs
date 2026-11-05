# T-Rex

The T-Rex model at the default [control parameters](details.md#control-parameters), as McStas
draws it: from above and along the beam. Both figures have the same scale on their two axes, in m.
[Technical details](details.md) lists every component, with its position, and describes each part.

## Seen from above

![T-Rex seen from above, z along the beam and x across it](figures/overview_top.png)

z is along the beam and x across it (projection `zx`). The circle at the origin is the drawing
of the source. The guides follow it in alternating colors up to `sample_pos`, 163.8 m along the
beam. Over that length the beam axis moves about 0.83 m in x, too little to see at this scale.
The short ticks across the beam are the choppers: `bw1` near 32 m, `bw2` near 40 m, the P-Chopper
near 108 m and the M-Chopper near 162 m.

This figure and the next need McStas, so they are not made with the site:
[`mcdgs.diagnostics.trex.overview`](../../reference/diagnostics.md) draws them with McStasScript's
`show_instrument` (mcdisplay), and `pixi run docs-mcstas` runs it. They show the instrument as of
that run.

## From detector to source

![T-Rex seen against the beam, x across it and y up](figures/overview_along_beam.png)

x is across the beam and y up (projection `xy`), with every component drawn on top of each other
and the view limited to x from −0.5 to 1.5 m and y from −1 to 1 m. The gray background is the
drawing of the source. The guides are the narrow band at y = 0, running from x ≈ 0 at the source
to x ≈ 0.83 m at the sample.

The chopper disks are 0.7 m across, each to one side of the beam. The two disks of `bw1` and
`bw2`, and of `ps1` and `ps2`, sit on the same side and overlap. `mc1` and `mc2` sit on
opposite sides:

| Disk | Choppers | Position | Side of the beam |
| ---- | -------- | -------- | ---------------- |
| Orange, left | `bw1`, `bw2` | x ≈ 0.01 to 0.03 m | Below |
| Red, middle | `ps1`, `ps2` | x ≈ 0.41 m | Above |
| Dark gray, right | `mc1` | x ≈ 0.82 m | Above |
| Light gray, right | `mc2` | x ≈ 0.82 m | Below |

The notches at the edge of the disks are the slits. The four of the P-Chopper are the clearest:
two of 20° and two of 35°. See [Chopper timing](details.md#chopper-timing) for how the mode puts
a narrow or a wide slit on the beam, and
[Simulate the choppers with tof](../../how-to/simulate-the-choppers-with-tof.md) for what the
choppers pass.
