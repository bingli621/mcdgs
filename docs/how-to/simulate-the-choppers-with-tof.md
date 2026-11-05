# Simulate the choppers with tof

[tof](https://scipp.github.io/tof) follows neutrons from the source through the choppers along a
single flight path. It knows nothing of guides or benders, so a million neutrons through T-Rex
take well under a second, with nothing to compile. It answers timing questions: which wavelengths
reach the sample, when they arrive, and how the mode or `rrm` changes that. For intensities, divergence or anything the
optics do, use [McStas](run-a-simulation.md).

The same instrument object feeds both. [`build_tof`](../reference/tof.md) turns it into a tof
model, with no change to the instrument.

## Build and run

```python
from mcdgs.instruments.trex.model import TRex
from mcdgs.tof.builder import build_tof

trex = TRex()
model = build_tof(trex, neutrons=1_000_000, seed=1)
result = model.run()
print(result)
```

The result says how many neutrons each chopper lets through and how many it stops:

```text
--8<-- "docs/instruments/trex/generated/tof_result.txt"
```

`result.plot()` draws the time-distance diagram, as on the [T-Rex](../instruments/trex/index.md#chopper-timing-in-tof)
page.

## What goes into the model

`build_tof` takes only what tof understands, by the type of the component:

| Component | In tof |
| --------- | ------ |
| `ESS_butterfly` | The source: tof's pulse for the beam port that `sector` and `beamline` name, from `Lmin` to `Lmax`, with `n_pulses` pulses at 14 Hz. W7 gets the T-Rex pulse and W3 the C-SPEC one; a beam port mcdgs does not know yet gets tof's standard ESS pulse. `cold_frac` and the focusing are not used. |
| `NXdisk_chopper` | A chopper, with the same slit edges, frequency, direction and phase as in McStas. |
| `Monitor_nD` | A detector, which records the neutrons reaching it. T-Rex has one, `sample_monitor`, at the sample position. With no monitor, the model has no detector. |
| Anything else, such as `Guide_gravity`, `Slit` or `Single_crystal` | Left out: tof takes the optics as letting everything through. |

A component's distance is the straight line from the moderator to its position. For T-Rex this
is within 1 mm of the path along the curved guide: 0.3 µs of flight at 2 Å, about 1 % of the
shortest chopper opening. The tof source itself is 5 cm from the moderator, where tof sampled
its pulse, so that its neutrons do not fly those 5 cm twice.

The choppers open and close exactly as in McStas. The tests check this every 0.5 µs over two
pulses, for every T-Rex disk in both modes, and check that a neutron reaches the sample in tof
exactly when the McStas disks would let it through. Two settings of `NXdisk_chopper` cannot be
carried over:

- A parked disk (`nu = 0`) raises an error: McStas holds it still at `park_angle`, while tof would
  take it as always open.
- The size of the openings, `jitter` and `abs_out` have no part in tof.

## Read the neutrons

Each detector holds every neutron of the source, with a mask for those stopped on the way:

```python
neutrons = result["sample_monitor"].data["pulse", 0]
reached = neutrons[~neutrons.masks["blocked_by_others"]]
print(reached.coords["wavelength"].min(), reached.coords["wavelength"].max())
print(reached.coords["toa"].to(unit="ms"))  # time of arrival
```

The coordinates also hold the `birth_time` and the `speed` of each neutron. The neutrons arrive
in sub-pulses, told apart by the gaps of about 4 ms between them.
[`summarize`](../reference/diagnostics.md) gives all of this at once:

```python
from mcdgs.diagnostics.trex.choppers import summarize

summary = summarize(result, "sample_monitor")
print(summary.neutrons, summary.sub_pulses, summary.first_arrival, summary.last_arrival)
```

## Compare settings

Building a model is cheap, so change the control parameters and build again:

```python
for mode in ("FF", "RR"):
    trex.update_params(mode=mode)
    reading = build_tof(trex, seed=1).run()["sample_monitor"]
    print(mode, reading)
```

With the same seed, both runs start from the same neutrons. At the default parameters, what
reaches the sample is:

--8<-- "docs/instruments/trex/generated/tof_modes.txt"

The wavelength band and the sub-pulses stay much the same, while the narrow slits of `"RR"` let
fewer neutrons through. The result above and this table are made from the current code for every
build of this site.
