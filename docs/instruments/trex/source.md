# Source

`TRexSource`, in
[`mcdgs.diagnostics.trex_01_source`](https://github.com/bingli621/mcdgs/blob/main/src/mcdgs/diagnostics/trex_01_source.py),
compares McStas and [tof](../../how-to/simulate-the-choppers-with-tof.md) at the entrance of
guide 1, the first optics of T-Rex. Nothing stands between the source and that point, so the two
should see the same pulse there: what differs is only how each models the ESS moderator.

## The setup

It is the smallest T-Rex: two components and no McStas instrument parameters.

| Component | What it is |
| --------- | ---------- |
| `source` | The T-Rex source, `ESS_butterfly` on beam port W7, one pulse at 2 MW, as in `TRex` but with the widest band, 0.1 to 10 Å, in place of the band around `lambda0`. |
| `source_monitor` | `source_monitor()` of [`mcdgs.instruments.trex.monitors`](../../reference/trex.md), a `Monitor_nD` that is the entrance of guide 1: its opening, 0.089955 m wide and 0.029992 m high, 1.895927 m from the source, 2.0 mm off its axis towards −x and turned by −0.215°, the place and rotation of `arm_nboa` in `TRex`. |

The monitor records one 2D map of time against wavelength:

| Axis | Range | Bins | Bin width |
| ---- | ----- | ---- | --------- |
| Time | 0.048 to 8.30 ms | 165 | 0.05 ms |
| Wavelength | 0.1 to 10 Å | 198 | 0.05 Å |

The wavelength axis is the band of the source. The time axis runs from the fastest neutron of
that band, 0.1 Å leaving at the start of the pulse, to the slowest, 10 Å leaving at its end, so it
holds every neutron of the source. Each axis has the fewest bins of its width that cover it. The
monitor has `restore_neutron=1`, so it does not change the beam.

`TRexSource` has no control parameters and no McStas parameters. It is defined like any instrument of mcdgs: a parameters class, then an `Instrument` subclass that
creates its components in beam order, with their positions.

<!-- fmt:off -->
```python
--8<-- "src/mcdgs/diagnostics/trex_01_source.py:trex_source"
```
<!-- fmt:on -->

## Running both models

The same `TRexSource` object runs in McStas and in tof. McStas writes the instrument file, runs it
and saves the monitor in its run folder, which McStasScript loads back:

<!-- fmt:off -->
```python
--8<-- "src/mcdgs/diagnostics/trex_01_source.py:run_mcstas"
```
<!-- fmt:on -->

tof builds its model from the same components: the `ESS_butterfly` becomes tof's pulse of beam
port W7 and the `Monitor_nD` a detector named `source_monitor`. It gives each neutron's time and
wavelength there. `read_comparison`, in
[`mcdgs.diagnostics.compare`](../../reference/diagnostics.md), reads the map of McStas in n/s per
bin and counts the tof neutrons on the same bins:

<!-- fmt:off -->
```python
--8<-- "src/mcdgs/diagnostics/trex_01_source.py:run_tof"
```
<!-- fmt:on -->

Both runs take the same `seed`. See [Run a McStas simulation](../../how-to/run-a-simulation.md)
and [Simulate the choppers with tof](../../how-to/simulate-the-choppers-with-tof.md) for each.

### Where the tof source sits

The two models place their source differently:

| | Source at | Monitor at | Flight to the monitor |
| --- | --- | --- | --- |
| McStas | 0 m, the moderator | 1.895927 m | 1.895927 m |
| tof | 0.05 m | 1.895927 m | 1.845927 m |

McStas emits each neutron at the moderator, with its time of emission. tof's pulse of beam port
W7, `ess-trex`, was not sampled there. It was sampled 5 cm away, so its times are those of the
neutrons when they pass 0.05 m, and they already include the flight over the first 5 cm. The
builder therefore keeps the tof source at 0.05 m, and the monitor at its distance from the
moderator, 1.895927 m. Both models then have the neutrons reach the monitor at the same time.

With the source at 0 m, tof would add the flight over those 5 cm a second time. That is
25 µs at 2 Å and 126 µs at 10 Å, half a bin to two and a half bins of 0.05 ms in the figure
below.

Every other distance in tof is the straight line from the moderator to the component. For
T-Rex that is within 1 mm of the path along the curved guide. Up to guide 1 the path is straight,
so here they are the same.

## Counts and intensity

Each simulated neutron stands for part of the real flux, its weight. The intensity of a bin is its
counts times the weight, in n/s. The counts depend only on how many neutrons are simulated; the
intensity is the physics, and should not.

- **McStas** sets the weights itself, so that the rays add up to the flux of the moderator. With
  ten times the `ncount`, a bin has ten times the rays, each with a tenth of the weight: the
  intensity stays the same, and only its noise goes down.
- **tof** draws its neutrons from a distribution with no flux: each has weight 1, so its
  intensity would be its counts, and would grow with the number of neutrons drawn. tof therefore
  takes its weight from McStas:

    ```text
    scale = McStas total (n/s) / tof total (counts)
    ```

    is the n/s that each tof neutron stands for. Scaled, tof has the level of McStas, and the shape
    in time and wavelength of its own model.

The comparison is thus of shapes, not of absolute flux: by construction both have the same total.

### What the scale depends on

The two totals behave differently with the number of neutrons. The McStas total is the flux Φ
through the opening of the monitor, in n/s, and does not depend on `ncount`. The tof total is the
number of tof neutrons that land in the bins: the `neutrons` drawn, times the fraction f of tof's
pulse that falls inside the window of the monitor. So

```text
scale ≈ Φ / (neutrons × f)
```

The scale is thus not a constant: with twice the tof neutrons it is half as large. What stays the
same is `scale × neutrons ≈ Φ / f`, the n/s of the whole tof pulse. That depends on:

- **The source model.** The power, 2 MW, the beam port, W7, the band, 0.1 to 10 Å, and how
  `ESS_butterfly` focuses its rays: Φ is proportional to the power, for example.
- **The acceptance of the monitor.** Its opening and its distance from the source set which of
  the rays of McStas pass, so which part of the flux Φ is. tof has no apertures: only its choppers
  remove neutrons, so here every tof neutron reaches the monitor, whatever its size. Move or
  resize the monitor and Φ, so the scale, changes.
- **The window of the bins.** The edges, 0.048 to 8.30 ms and 0.1 to 10 Å, set both which part of Φ is
  counted and f. If tof's pulse reaches beyond 10 Å, f is below 1. At the defaults of `main`,
  8.46e11 n/s and 8.46e5 n/s per tof neutron make 1.00e6 neutrons in the bins, as many as were
  drawn: here f is 1, to three figures.
- **The statistics.** Both totals are noisy, so the scale has a relative error of about
  √(1/(neutrons × f) + (σ/Φ)²), with σ the error of the McStas total. That is 0.1 % at a million
  tof neutrons, but 1 % at 10,000, and it is not in the error bars of tof.

Since tof has no flux, the scale is set where both models see the same beam, here at the
source. Further down, keep that same scale rather than computing a new one at each monitor:
scaled again to the total of McStas there, tof would agree with it by construction, and a
difference in what the guides and choppers let through would not show. Keep it with the same
number of tof neutrons; with another number, keep `scale × neutrons`, as the scale goes as
1/neutrons. `read_comparison` takes it as `scale`; left out, it scales to the McStas total at
that monitor.

## Error bars

A bin of N neutrons has the Poisson error √N, a relative error of 1/√N. In n/s it is √N times the
weight. With a hundred times the neutrons, the error in counts is ten times larger, but the
relative error is a tenth, and so is the error in n/s, as the weight is a hundredth. The error
bars shrink with more neutrons in both models, while the intensity stays the same.

McStas gives the error of each bin itself. For tof it is √N of the counts, scaled by the same
factor as the counts. The bins are independent, so when the maps are summed over one axis their
errors add in quadrature.

## The figure

`main` runs both models with the same seed, reads each monitor on the same bins, prints the mean
time and wavelength of each and the McStas total, scales tof to McStas, and saves
`trex_source.png` next to the `mccode.h5` of the McStas run, or to the folder given with `--out`.
The figure has two rows:

- **Top:** the time and wavelength map of McStas and of tof, side by side on one color scale, in
  n/s per 0.05 ms and 0.05 Å.
- **Bottom:** the wavelength spectrum, summed over time, and the time profile, summed over
  wavelength, of both models on the same axes, with their error bars. The legend gives the
  McStas `ncount`, the number of tof neutrons and the n/s that each of them stands for.

![McStas and tof at the entrance of guide 1: time and wavelength maps, the wavelength spectrum and the time profile](figures/trex_source.png)

This is `main` at its defaults: 1e6 McStas neutrons, 1e6 tof neutrons and seed 1. The total at
the entrance of guide 1 is 8.46e11 n/s, so each tof neutron stands for 8.46e5 n/s. Both models
have the same means, 2.84 ms and 2.40 Å.

The two models agree within their error bars almost everywhere. The spectrum has the thermal
peak near 1.1 Å, the cold one near 2.6 Å, and a third peak at the shortest wavelengths, below
0.3 Å. Only there do they differ: in the first bins McStas is about 15 % higher than tof. In time,
both rise over the 2.86 ms of the pulse and peak near 2.8 ms; the long wavelengths arrive later,
up to about 8 ms at 10 Å.

The figure needs McStas, so it is not made with the site: `pixi run docs-diagnostics` makes it
again, with the other files of the T-Rex pages.

To make it yourself, run:

```bash
pixi run python -m mcdgs.diagnostics.trex_01_source
```

Call `main` with other numbers for a finer figure:

```bash
pixi run python -c "from mcdgs.diagnostics.trex_01_source import main; main(ncount=1e8, neutrons=10_000_000)"
```

## Estimating the count rate of a monitor

The intensity on this monitor is also what a real beam monitor at that place would see, so the same
map tells you which monitor can stand there. What decides it is the rate the monitor counts,
above all at the peak of each pulse, not the average flux.

### Per bin, per second, per pulse

A value of the map is in n/s per bin: the neutrons per second that pass the opening of the
monitor, with a time of arrival in that 0.05 ms and a wavelength in that 0.05 Å. The time is the
time within each pulse; the second is the average over the 14 Hz pulses. Summed over all bins, it
is the total through the opening:

| Quantity | How | At the entrance of guide 1 |
| -------- | --- | -------------------------- |
| Total rate | sum of all bins | 8.46e11 n/s |
| Per pulse | total / 14 Hz | 6.0e10 n per pulse |
| Per area | total / (8.9955 cm × 2.9992 cm) | 3.1e10 n/s/cm² |

Halving the widths of the bins halves each value of the map, but leaves these sums unchanged;
only a change of the window changes them. Divide a bin by its widths, 0.05 ms and 0.05 Å, for a
density in n/s/ms/Å that does not depend on the bins.

### The peak rate

Each pulse reaches the monitor within a few ms of the 71.4 ms between pulses, so while it passes
the rate is far above the average. From the time profile, the sum of the map over wavelength, the
instantaneous rate in a time bin is

```text
instantaneous rate (n/s) = time profile (n/s per bin) / 14 Hz / 0.05 ms
```

Its maximum is the peak rate. At the entrance of guide 1, in the map of `main` at its defaults,
it is 2.0e13 n/s at 2.7 ms, some 24 times the average. That is what the length of the pulse
suggests: most of the 6.0e10 neutrons of a pulse arrive within about 3 ms, and
6.0e10 / 3 ms ≈ 2e13 n/s.

### What the monitor counts

A beam monitor must leave the beam almost unchanged, so it detects only a small fraction ε of the
neutrons, typically 1e-5 to 1e-3. For ³He, ¹⁰B and most gas monitors ε grows about as the
wavelength, so a cold beam counts more than its flux alone suggests. With ε(λ) from the data sheet,
on the centers of the wavelength bins:

```python
time_profile = (mcstas_intensity * efficiency).sum(axis=1)  # counts/s per time bin
average_rate = time_profile.sum()  # counts/s
peak_rate = time_profile.max() / 14 / 0.05e-3  # counts/s, during the pulse
```

`mcstas_intensity` is the map of `main`, of shape (time bins, wavelength bins), and `efficiency`
is ε at the wavelength centers.

### Choosing

Compare the peak rate with what the monitor can count. While peak rate × τ is small, a counting
monitor with a dead time τ loses a fraction of about peak rate × τ of its counts; keep it below
0.1, or the monitor no longer counts in proportion to the beam. Near 1 and above, it saturates:
it cannot count much faster than 1/τ, whatever the beam. If peak rate × τ is too large, take a
lower ε, a smaller opening, or a faster monitor. Leave a factor of 2 to 10 for a higher source
power and for the model.

For example, with ε = 1e-4 at 1.8 Å, growing as λ: at the mean wavelength of 2.4 Å ε is 1.3e-4,
so the monitor counts about 1.1e8 counts/s on average and, from the map, 2.5e9 counts/s at the
peak. With a dead time of 100 ns, peak rate × τ is about 250, far above 0.1: the monitor is
saturated, at most 1e7 counts/s. No counting monitor of that ε works at the entrance of guide 1,
which is close to the moderator and inside the monolith. Monitors stand further downstream, after
the choppers, where the beam is weaker and shorter in wavelength band. Put the `Monitor_nD`
there, with the opening of the real monitor, and repeat.

McStas leaves out the fast neutrons and the γ rays from the source, which a monitor near the
source also sees. They do not change the rate of slow neutrons, but they decide which kinds of
monitor can stand there.
