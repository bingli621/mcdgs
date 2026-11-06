# Extraction guides

`Extraction`, in
[`mcdgs.diagnostics.trex.extraction`](https://github.com/bingli621/mcdgs/blob/main/src/mcdgs/diagnostics/trex/extraction.py),
is T-Rex up to the beam monitor after guide 4: the source, the
[extraction guides](details.md#extraction-guides) 1 to 3, and guide 4 with its
[two paths](details.md#two-paths-through-guide-4), the thermal guide 4a with the bender out and
the bender and the cold guide 4b with it in. It has the same components as `TRex`, with the same
values, but for the band: its source keeps the widest band, 0.1 to 10 Å, and its beam monitor
that band too, while those of `TRex` follow `lambda0`.

`main` runs McStas once per path and tof once, and makes two figures of each path: McStasScript's
diagnostics along the guides, and McStas and tof side by side at the beam monitor. The figures are
`main` at its defaults: 1e6 McStas neutrons per path, 1e6 tof neutrons and seed 1.

## Beam along the guides

Each row is one position in the beam: after the source, after guide 3 and at the exit of guide 4b.
The columns are the beam profile, wavelength against x, and the spectrum in n/s per bin. x is
across the beam in the frame of the component, so after guide 3 it spans the exit of guide 3,
±31.9 mm, and after guide 4b that of guide 4, ±30.1 mm.

**Bender out**, through the thermal guide 4a:

![Beam after the source, guide 3 and guide 4 with the bender out](figures/trex_extraction_bender_out_diagnose.png)

**Bender in**, through the bender and the cold guide 4b:

![Beam after the source, guide 3 and guide 4 with the bender in](figures/trex_extraction_bender_in_diagnose.png)

The top two rows are the same in both figures. After the source the two moderators are side by
side: the thermal one at x < 0 with its peak near 1.2 Å, the cold one at x > 0 with its peak near
2.6 Å. After guide 3 the two are mixed across the width of the guide, and the spectrum still has
both peaks. The bottom rows differ: without the bender both peaks reach the exit, while with it
the thermal peak is mostly gone and the cold peak stays.

The maps have 100 × 100 bins, so at 1e6 neutrons each bin holds only a few rays and the maps are
grainy. The spectra, on one axis, are smooth.

## McStas and tof at the beam monitor

The beam monitor is in the gap after guide 4, at 5.902 m, 0.06 m wide and 0.0463 m high. It
records one map of time, 0.15 to 18.4 ms in 366 bins of 0.05 ms, against wavelength, 0.1 to 10 Å
in 495 bins of 0.02 Å. The time runs from the fastest neutron of the band, 0.1 Å, leaving at the
start of the pulse, to the slowest, 10 Å, leaving 3.5 ms later. Both models see the same monitor;
see [Source](source.md) for how the two are read and compared. With 1e6 neutrons spread over
that many bins, the McStas map and spectra are grainy.

tof has no guides and no bender: every neutron of its pulse reaches the monitor, at the time of
its straight flight from the source. Its counts are scaled to the McStas total at the monitor, on
each path, so both have the same total and the figure compares their shapes. What differs is
what guides 1 to 4 do to the beam.

**Bender out**, 1.48e11 n/s at the monitor, 1.48e5 n/s per tof neutron:

![McStas and tof at the beam monitor with the bender out](figures/trex_extraction_bender_out.png)

**Bender in**, 1.04e11 n/s at the monitor, 1.04e5 n/s per tof neutron:

![McStas and tof at the beam monitor with the bender in](figures/trex_extraction_bender_in.png)

With the bender out, McStas keeps both peaks of the source but has fewer short wavelengths and
more long ones than tof: the guides pass a long wavelength at a larger angle, so more of it
reaches the monitor. Below 0.3 Å McStas has a fraction of tof's peak. With the bender in, McStas
has almost nothing below about 1.5 Å and peaks near 2.5 Å, and its time profile peaks about
1.5 ms later than tof's, since the beam is made of longer, slower wavelengths.

## Making the figures

The figures need McStas, so they are not made with the site: `pixi run docs-diagnostics` makes
them again, with the other files of the T-Rex pages. To make them yourself, run:

```bash
pixi run python -m mcdgs.diagnostics.trex.extraction
```

Each path's two figures go to its McStas run folder, `mcstas_output/trex_extraction_bender_out/`
and `..._in/`. That folder also holds McStasScript's event lists of the three diagnostics
positions, about 190 MB for each path at 1e6 neutrons. Call `main` with other numbers for finer
figures:

```bash
pixi run python -c "from mcdgs.diagnostics.trex.extraction import main; main(ncount=1e7)"
```
