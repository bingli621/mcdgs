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
ferromagnet with exchange $J$ between nearest neighbors, spins $S$ and an anisotropy gap $\Delta$:

$$
\hbar\omega(\mathbf{q}) = \Delta + 4JS\,(3 - \cos q_x a - \cos q_y a - \cos q_z a)
$$

With $a$ = 4 Å, $\Delta$ = 5 meV and $JS$ = 0.625 meV, the band runs from the 5 meV gap at the zone
center to 20 meV at the zone corner, $(\pi/a)(1, 1, 1)$. A powder sees, at each |Q|, every
wavevector of that length. `write_magnon_sqw` averages the magnon over 20,000 directions on a
Fibonacci lattice of the sphere and smooths each row with a Gaussian of $\sigma$ = 0.03 meV, far below
the resolution of T-Rex. At small |Q| the powder band is a narrow dispersive line, $\Delta + 2JS\,(\lvert Q \rvert a)^2$.
At larger |Q| it is a band between a lowest and a highest energy, `magnon_band`, with sharp
edges.

| | |
| - | - |
| Sample | a solid cylinder 20 mm across and 10 mm high, the size of the vanadium, `Isotropic_Sqw` |
| S(\|Q\|, E) | 0 to 10 Å⁻¹ in steps of 0.02 Å⁻¹, −20.2 to 20.2 meV in steps of 0.02 meV, symmetric in E, each side of each row normalized to one magnon |
| Cross-sections | 5 b coherent per 64 Å³ cell, only the scale of the intensity; no absorption, no incoherent scattering |
| Taken as it is | `T=-1`: no detailed balance; `classical=0`, `norm=1`; the table holds both sides, so the magnon is annihilated as often as it is created, a Bose factor of 1 |
| Scattering | single scattering only (`order=1`), within `d_phi` = 52° of the horizontal plane, 90% of the neutrons made to scatter (`p_interact=0.9`), each neutron sent through it 4 times (`SPLIT 4`) |
| Detector | the banana detector, as for vanadium |

The magnetic form factor and the polarization factor are left out, so the intensity does not
fall with |Q|. `Isotropic_Sqw` has no horizontal focusing: it scatters over the whole horizontal
circle, and the banana detector catches the half of it between −36° and 144°.

### Both sides of E, with a Bose factor of 1

A magnon is created when the neutron loses energy, at positive E, and annihilated when it gains
energy, at negative E. In a sample at a temperature $T$, detailed balance ties the two sides by
the Bose factor:

$$
S(Q, -E) = e^{-E / k_B T}\, S(Q, E)
$$

Here the two sides are equal at every energy transfer, a Bose factor of 1, the limit
$k_B T \gg E$. `write_magnon_sqw` writes the table on both sides itself, $S(Q, -E) = S(Q, E)$,
and `T=-1` makes `Isotropic_Sqw` take it as it is, with no detailed balance. The figure then shows
the band on both sides of E, as strong on each: any difference between them comes from the
instrument or the reduction, not from the sample.

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
`sample_pos` with 5e9 rays, then from that beam through the powder to the detector. Each neutron of
the beam goes through the powder 4 times, which gives 67 million detector events.

`reduce_runs` is the reduction of the vanadium: the sub-pulses of `subpulse_timing`, the pixels
past 5° from the beam, the neutrons of each sub-pulse divided by the solid angle of their pixel
and by its flux, and `q_e_maps`, S(Q, E) normalized by the coverage of each bin, in bins of
0.25 meV and 0.1 Å⁻¹.

The intensity a detector measures is the double-differential cross section, which holds a factor
$k_f/k_i$ on top of S(Q, E):

$$
\frac{d^2\sigma}{d\Omega\, dE_f} = \frac{\sigma}{4\pi}\, \frac{k_f}{k_i}\, S(Q, E),
\qquad \frac{k_f}{k_i} = \sqrt{\frac{E_f}{E_i}} = \sqrt{1 - \frac{E}{E_i}}
$$

`Isotropic_Sqw` scatters with it, so a neutron that gains energy, and leaves faster, counts more
than one that loses as much. `q_e_maps` removes it: it divides each neutron by its $k_f/k_i$,
from `kf_over_ki`, before it histograms them, so the maps are S(Q, E) itself. Without that, the
gain side of a sub-pulse would be brighter than its loss side, though the sample scatters the same
on both: by about 2 across the band at the lowest Ei, and by up to 3 at the edges of an arrival
window, $\sqrt{1.8 / 0.2}$ at $E = \mp 0.8\, E_i$.

## S(Q, E)

`plot_magnon_maps` draws the 15 maps, from the highest Ei to the lowest, from −2 to 23 meV, with
`plot_q_e_maps` of the vanadium: one linear color scale for all of them, white where the
sub-pulse cannot reach. Over each map, dashed white lines show the lowest and the highest energy
of the magnon at each |Q|, without resolution. How far the intensity spreads past them is the
resolution of that sub-pulse.

![S(Q, E) of the powder magnon in each of the 15 sub-pulses in FF, |Q| across and E up, from the highest Ei to the lowest, normalized by the coverage of each bin, with the edges of the band of the model as dashed white lines](figures/magnon-powder/trex_magnon_powder_ff_sqe.png)

The powder magnon is not a smooth band. At every zone center, $\lvert Q \rvert = (2\pi/a)\sqrt{h^2 + k^2 + l^2}$ =
1.57, 2.22, 2.72, 3.14 Å⁻¹, ..., the lower edge drops back to the gap, and a dark V of little
intensity rises from it. Between the Vs, bright arches run across the band. At small |Q| the
steep dispersion, from 5 meV to about 10 meV below 1 Å⁻¹, is the left edge of the intensity.
How sharply each sub-pulse draws this pattern is its resolution:

- **The highest Ei, 47 to 36 meV (0 to 2), see the whole band but blur it.** The Vs and the
  arches merge into one broad blob, from about 8 to 17 meV, and the gap at 5 meV is a soft
  floor. Sub-pulse 0, in the tail of the source spectrum, has the fewest neutrons: it is the
  noisiest, and its hottest bins set the top of the color scale.
- **From 23 to 14 meV (6 to 11), the pattern is sharp.** The Vs, the arches, the gap and the
  dispersion at small |Q| all have clear edges, though these sub-pulses see only the band below
  0.8 Ei.
- **The lowest Ei, 13 to 11 meV (12 to 14), see only the bottom of the band**, up to about
  9 meV for the 11.4 meV sub-pulse, and are faint and grainy. About a fifth of the neutrons
  that should scatter are transmitted instead, as McStas reports, because it finds no (|Q|, E)
  of the table that sends them towards the detector. That loss is heaviest where little of the
  band can be reached, so the intensity of these sub-pulses is somewhat too low.

Below the gap, from about 9 to 14, there is faint intensity that the model does not have. It is
frame overlap: neutrons of one sub-pulse that lose enough energy in the powder reach the detector
in the arrival window of the next sub-pulse, so the reduction gives them the Ei of that one, and
so a wrong energy transfer.

## Making the figure

```console
pixi run python -c "from mcdgs.diagnostics.trex.magnon_powder import main; main()"
```

This runs McStas, then `analyze`, which writes `trex_magnon_powder_ff_sqe.png` into the results
folder. `pixi run python -m mcdgs.diagnostics.trex.magnon_powder` only reduces the runs again,
without McStas.
