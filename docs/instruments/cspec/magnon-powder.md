# Powder magnon

`CSPECMagnonPowder`, in
[`mcdgs.diagnostics.cspec.magnon_powder`](https://github.com/bingli621/mcdgs/blob/main/src/mcdgs/diagnostics/cspec/magnon_powder.py),
replaces the vanadium of [Elastic resolution](elastic-resolution-incoherent-short.md) with a powder
that scatters from one magnon branch. It draws S(Q, E) for each of the 14 sub-pulses in one
figure. Vanadium gives one line at E = 0, so it measures the resolution at a single energy
transfer. The magnon spreads intensity across E and |Q|, and its band has sharp edges, so the
figure shows how each sub-pulse blurs a dispersive excitation over its whole Q–E window. It
follows [Powder magnon of T-Rex](../trex/magnon-powder.md), with a magnon scaled to the incident
energies of CSPEC.

The reduction and the figure are those of T-Rex. The parts both instruments share are in
[`mcdgs.diagnostics.magnon_powder`](../../reference/diagnostics.md).

## The sample

The powder is `magnon_powder`, an `Isotropic_Sqw`, in `samples/inelastic/powder.py`, as on T-Rex.
It reads a table of S(|Q|, E) that `write_magnon_sqw` writes. The magnon is that of a
simple-cubic ferromagnet with exchange $J$ between nearest neighbors, spins $S$ and an anisotropy
gap $\Delta$:

$$
\hbar\omega(\mathbf{q}) = \Delta + 4JS\,(3 - \cos q_x a - \cos q_y a - \cos q_z a)
$$

The sub-pulses of CSPEC have Ei from 15.1 to 5.7 meV, and their arrival windows reach an energy
transfer of 0.8 Ei, from 12 to 4.6 meV. The magnon of T-Rex, from 5 to 20 meV, is mostly out of
that reach, so CSPEC has its own. It keeps the cell of T-Rex, $a$ = 4 Å, so the zone centers are
at the same |Q|. With $\Delta$ = 1 meV and $JS$ = 7/24 meV = 0.292 meV, the band runs from the
1 meV gap at the zone center to 8 meV at the zone corner, $(\pi/a)(1, 1, 1)$. A powder sees, at
each |Q|, every wavevector of that length. `write_magnon_sqw` averages the magnon over 20,000
directions on a Fibonacci lattice of the sphere and smooths each row with a Gaussian of
$\sigma$ = 0.01 meV, a quarter of the narrowest FWHM of CSPEC, 0.10 meV, from
[Elastic resolution](elastic-resolution-incoherent-short.md). At small |Q| the powder band is a narrow
dispersive line, $\Delta + 2JS\,(\lvert Q \rvert a)^2$. At larger |Q| it is a band between a
lowest and a highest energy, `magnon_band`, with sharp edges.

| | |
| - | - |
| Sample | a solid cylinder 20 mm across and 10 mm high, the size of the vanadium, `Isotropic_Sqw` |
| S(\|Q\|, E) | 0 to 6 Å⁻¹ in steps of 0.02 Å⁻¹, −8.06 to 8.06 meV in steps of 0.005 meV, symmetric in E, each side of each row normalized to one magnon |
| Cross-sections | 5 b coherent per 64 Å³ cell, only the scale of the intensity; no absorption, no incoherent scattering |
| Taken as it is | `T=-1`: no detailed balance; `classical=0`, `norm=1`; the table holds both sides, so the magnon is annihilated as often as it is created, a Bose factor of 1 |
| Scattering | single scattering only (`order=1`), within `d_phi` = 55.1° of the horizontal plane, the height of the window of the vanadium, 90% of the neutrons made to scatter (`p_interact=0.9`), each neutron sent through it 4 times (`SPLIT 4`) |
| Detector | the twelve detector panels, as for vanadium |

The magnetic form factor and the polarization factor are left out, so the intensity does not
fall with |Q|. `Isotropic_Sqw` has no horizontal focusing: it scatters over the whole horizontal
circle, and the panels catch the part of it between −31.5° and 143.8°.

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

`simulate` writes the table into the builds folder of the workspace, `sqw_path`, as
`cspec_magnon_powder.sqw`, apart from the table of T-Rex, before it runs McStas. McStas reads it
when the instrument starts, in both parts of the run, so the path is absolute.

<!-- fmt:off -->
```python
--8<-- "src/mcdgs/diagnostics/cspec/magnon_powder.py:instrument"
```
<!-- fmt:on -->

## Runs and reduction

`simulate` runs the default control parameters, `lambda0` = 3 Å, `nu_mc` = 224 Hz, with the
masked source, in the two parts of
[Elastic resolution](elastic-resolution-incoherent-short.md#two-runs-split-at-the-sample):

1. **From the source to `sample_pos`**, with 5e9 rays, in 17.5 minutes.
2. **From that beam through the powder to the panels**, in 22 seconds. Each neutron of the beam
   goes through the powder 4 times, which gives 3.9 million detector events.

`reduce_runs` is the reduction of the vanadium: the sub-pulses of `subpulse_timing` at
`mc_monitor`, every pixel of the panels, the neutrons of each sub-pulse divided by the solid
angle of their pixel and by its flux, and `q_e_maps`, S(Q, E) normalized by the coverage of each
bin, in bins of 0.1 meV and 0.05 Å⁻¹.

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

`plot_magnon_maps` draws the 14 maps, from the highest Ei to the lowest, over the whole arrival
windows: from −13 to 13 meV, $\mp 0.8$ of the highest Ei, so both the energy the neutron gains and
the energy it loses. It draws them with `plot_q_e_maps` of the vanadium: one linear color scale
for all of them, white where the sub-pulse cannot reach. The scale reaches the 98th percentile of
the bins, lower than the 99.9th of T-Rex, so that the band shows in the upper colors rather than
only the hottest bins of sub-pulse 13. Over each map, dashed white lines show the lowest and the
highest energy of the magnon at each |Q|, without resolution, at positive E and mirrored at
negative E. How far the intensity spreads past them is the resolution of that sub-pulse. The
sample puts nothing between the two bands nor beyond them, so whatever lies there comes from the
instrument or the reduction.

![S(Q, E) of the powder magnon in each of the 14 sub-pulses, |Q| across and E up, from the highest Ei to the lowest, over the energy gain and the energy loss, normalized by the coverage of each bin, with the edges of the band of the model as dashed white lines on both sides](figures/magnon-powder/cspec_magnon_powder_sqe.png)

The pattern is that of T-Rex, at the same |Q| and a smaller energy scale. At every zone center,
$\lvert Q \rvert = (2\pi/a)\sqrt{h^2 + k^2 + l^2}$ = 1.57, 2.22, 2.72, 3.14 Å⁻¹, ..., the inner
edge drops back to the gap, and a dark V of little intensity rises from it. Between the Vs,
bright arches run across the band. At small |Q| the steep dispersion, from 1 meV to about 3 meV
at 0.8 Å⁻¹, is the bright inner edge of the intensity: there the whole band is one narrow line,
so its intensity per meV is the highest of the map. The same pattern stands at negative E,
upside down.

- **Every sub-pulse draws the pattern sharply.** The FWHM of CSPEC, from 0.6 meV at 15 meV to
  0.1 meV at 5.7 meV, is small next to the 7 meV band, so even the highest Ei, 15 to 12.5 meV
  (0 to 2), keep the Vs, the arches and the gap apart. Their edges are the softest, and they
  sharpen as Ei falls.
- **The highest Ei, 15 to 12.5 meV (0 to 2), see the whole band on both sides**, up to |Q| of
  about 4 Å⁻¹. The lower the Ei, the less of it: sub-pulse 6, at 9.1 meV, reaches 7.3 meV on
  each side, and sub-pulse 12, at 6.1 meV, 4.9 meV and 3 Å⁻¹.

The two sides test the instrument against each other, since the sample scatters the same on
both. Inside the band, over the (|Q|, |E|) both sides reach, the mean S(Q, E) of the gain side is
0.96 to 0.99 of that of the loss side for sub-pulses 0 to 7. They part where the instrument adds
what the sample does not have, more than on T-Rex: the band reaches 8 meV, most of the Ei of the
lower sub-pulses, which follow each other every 4.5 ms.

- **Frame overlap from earlier sub-pulses, from sub-pulse 7 on.** A neutron that loses 5 to 8 meV
  is slow enough to reach the detector in the arrival window of a later sub-pulse, so the
  reduction gives it the Ei of that one, and so a wrong energy transfer, anywhere across the
  window. Between the two bands, at least one FWHM of the elastic line from their edges, there is
  at most 0.1 % of the intensity of the band up to sub-pulse 6, then 0.7 % to 2.4 % for 7 to 11
  and 6 % for 12. These neutrons add to the loss side, and the gain side falls behind: 0.86, 0.74,
  0.47 and 0.50 of it for 8 to 11.
- **Frame overlap from later sub-pulses, at the top of the window.** A neutron of the next
  sub-pulse that gains 3 to 8 meV is fast enough to reach the detector at the very end of the
  window of the one before. The reduction gives it an energy loss of about 0.8 Ei, where a
  millisecond of arrival time spans the fewest meV, so these few neutrons pile up in the top
  bins. That is the bright stripe along the top of the maps of 9 to 11, with 6 to 11 times the
  intensity of the bins just below it; in 8 the stripe is 2 times, from late neutrons of 7 alone.
- **Sub-pulse 13, the last at 5.7 meV, is mostly frame overlap.** Between the bands it has 32 %
  of the intensity of its band, and its gain side only 0.52 of its loss side. Its arrival window
  overlaps that of sub-pulse 12, and a neutron in both goes to 13, so it also takes the neutrons
  of 12 that lose the most energy.
- **Intensity past the gain band.** Beyond the outer edge of the band, at least one FWHM of the
  elastic line from it, the loss side holds at most 0.7 % of the band up to sub-pulse 5, the gain
  side 1 % to 4 % away from the edges of the coverage. It is most likely the resolution, which
  broadens with Ef, so the gain band is blurred more than the loss band; its width has not been
  measured yet. At |Q| below 1.2 Å⁻¹ the gain side holds more, 5 % to 20 % of the band, in a
  bright streak past the band along the low-angle edge of the coverage, where each bin sees little
  of the detector. It has not been explained yet.

Fewer of the neutrons that should scatter are lost than on T-Rex: McStas reports that 3.4 % of
them, 312,736 of 9.1 million, are transmitted instead, because it finds no (|Q|, E) of the table
that sends them towards the detector.

## Making the figure

```console
pixi run python -c "from mcdgs.diagnostics.cspec.magnon_powder import main; main()"
```

This runs McStas, then `analyze`, which writes `cspec_magnon_powder_sqe.png` into the results
folder. `pixi run python -m mcdgs.diagnostics.cspec.magnon_powder` only reduces the runs again,
without McStas.
