# Powder magnon

`TRexMagnonPowder`, in
[`mcdgs.diagnostics.trex.magnon_powder`](https://github.com/bingli621/mcdgs/blob/main/src/mcdgs/diagnostics/trex/magnon_powder.py),
replaces the vanadium of [Elastic resolution](elastic-resolution-incoherent.md) with a powder
that scatters from one magnon branch. It draws S(Q, E) for each of the 15 sub-pulses in one
figure. Vanadium gives one line at E = 0, so it measures the resolution at a single energy
transfer. The magnon spreads intensity across E and |Q|, and its band has sharp edges, so the
figure shows how each sub-pulse blurs a dispersive excitation over its whole Q–E window.

## The sample

The powder is `magnon_powder`, an `Isotropic_Sqw`, in `samples/inelastic/powder.py`. It reads a
table of S(|Q|, E) that `write_magnon_sqw` writes. The magnon is that of a simple-cubic
ferromagnet with exchange J between nearest neighbors, spins S and an anisotropy gap:

ħω(**q**) = Δ + 4JS (3 − cos qₓa − cos q_y a − cos q_z a)

With a = 4 Å, Δ = 5 meV and JS = 0.625 meV, the band runs from the 5 meV gap at the zone
center to 20 meV at the zone corner, (π/a)(1, 1, 1). A powder sees, at each |Q|, every
wavevector of that length. `write_magnon_sqw` averages the magnon over 20,000 directions on a
Fibonacci lattice of the sphere and smooths each row with a Gaussian of σ = 0.03 meV, far below
the resolution of T-Rex. At small |Q| the powder band is a narrow dispersive line, Δ + 2JS(|Q|a)².
At larger |Q| it is a band between a lowest and a highest energy, `magnon_band`, with sharp
edges.

| | |
| - | - |
| Sample | a solid cylinder 20 mm across and 10 mm high, the size of the vanadium, `Isotropic_Sqw` |
| S(\|Q\|, E) | 0 to 10 Å⁻¹ in steps of 0.02 Å⁻¹, 0 to 20.2 meV in steps of 0.02 meV, each row normalized to one magnon |
| Cross-sections | 5 b coherent per 64 Å³ cell, only the scale of the intensity; no absorption, no incoherent scattering |
| Taken as it is | `T=-1`: no detailed balance, the magnon is only created; `classical=0`, `norm=1` |
| Scattering | single scattering only (`order=1`), within `d_phi` = 52° of the horizontal plane, half of the neutrons made to scatter (`p_interact=0.5`) |
| Detector | the banana detector, as for vanadium |

The magnetic form factor and the polarization factor are left out, so the intensity does not
fall with |Q|. `Isotropic_Sqw` has no horizontal focusing: it scatters over the whole horizontal
circle, and the banana detector catches the half of it between −36° and 144°.

`simulate` writes the table into the builds folder of the workspace, `sqw_path`, before it runs
McStas. McStas reads it when the instrument starts, in both parts of the run, so the path is
absolute.

<!-- fmt:off -->
```python
--8<-- "src/mcdgs/diagnostics/trex/magnon_powder.py:instrument"
```
<!-- fmt:on -->

## Runs and reduction

`simulate` runs FF only, the default `lambda0` = 2 Å and `rrm` = 18, with the masked source, in
the two parts of [Diffraction](diffraction.md#two-runs-split-at-the-sample): from the source to
`sample_pos` with 1e9 rays, then from that beam through the powder to the detector.

`reduce_runs` is the reduction of the vanadium: the sub-pulses of `subpulse_timing`, the pixels
past 5° from the beam, the neutrons of each sub-pulse divided by the solid angle of their pixel
and by its flux, and `q_e_maps`, S(Q, E) normalized by the coverage of each bin, in bins of
0.1 meV and 0.05 Å⁻¹.

## S(Q, E)

`plot_magnon_maps` draws the 15 maps, from the highest Ei to the lowest, from −2 to 23 meV, with
`plot_q_e_maps` of the vanadium: one linear color scale for all of them, white where the
sub-pulse cannot reach. Over each map, dashed white lines show the lowest and the highest energy
of the magnon at each |Q|, without resolution. How far the intensity spreads past them is the
resolution of that sub-pulse.

- **The highest Ei see the whole band**, 5 to 20 meV, and blur it the most: their resolution
  is several meV.
- **The lowest Ei see only its bottom.** Their arrival window ends at 0.8 Ei, about 9 meV for
  the 11.4 meV sub-pulse, but they resolve the gap and the dispersion at small |Q| best.

## Making the figure

```console
pixi run python -c "from mcdgs.diagnostics.trex.magnon_powder import main; main()"
```

This runs McStas, then `analyze`, which writes `trex_magnon_powder_ff_sqe.png` into the results
folder. `pixi run python -m mcdgs.diagnostics.trex.magnon_powder` only reduces the runs again,
without McStas.
