# Read the output

Every run writes a McStas NeXus (HDF5) file, `mccode.h5`. A scan has one entry per value, `entry1`,
`entry2` and so on, and a single run has one. In the examples `run` and `scan` are what
`simulation.run` and `simulation.scan` returned, or what
[`latest_run`](run-a-simulation.md#find-and-clean-up-runs) found: `data_file` is the file and
`folder` the run folder.

## A monitor

Monitors are read with `h5py`. This reads `sample_monitor`, the `Monitor_nD` of T-Rex with
`filename="sample_monitor"`, as in [Install and first run](../getting-started.md):

```python
import h5py

with h5py.File(run.data_file) as file:
    components = file["entry1/instrument/components"]
    monitor = next(
        c for name, c in components.items() if name.partition("_")[2] == "sample_monitor"
    )
    time_map = monitor["output/sample_monitor_t_L"]
    times = time_map["TOF__s_"][()]  # bin centers, in s
    wavelengths = time_map["Wavelength__Angs_"][()]  # bin centers, in Å
    intensity = time_map["data"][()]  # one row per time bin, one column per wavelength bin
    error = time_map["errors"][()]
```

Component groups are numbered in beam order, such as `0076_sample_monitor`, which is why this
finds the group by the name after the number: matching its end would also find
`0073_pre_sample_monitor`. McStas adds a suffix for the variables to the file name, `_t_L` for
time and wavelength. To see what a monitor wrote, use `list(monitor["output"])`.

Next to the `mccode.h5` of a scan, `mcrun` writes `mccode.dat`: for each value, the total
intensity and error of every histogram monitor.

To read every entry of a scan, count through the entries, since sorted by name `entry10` comes
before `entry2`:

```python
with h5py.File(scan.data_file) as file:
    for n in range(1, len(file) + 1):
        entry = file[f"entry{n}"]
```

!!! warning "McStasScript reads only the first entry of a scan"

    `ms.load_data(scan.folder)` reads `entry1` alone, the first value, without saying so. It
    suits a single run; read a scan by its entries, as above, or with `read_mcstas_map` below.

### Into scipp

[`read_mcstas_map`](../reference/mcstas.md) reads the map of a monitor of a run into
[scipp](https://scipp.github.io), in n/s per bin, with the square of its errors as variances, on
the bin edges of the monitor, which [`monitor_bins`](../reference/reduction.md) gives from its
`bin_limits`. It takes the record of the run: a single run it reads as McStasScript loads it, a
scan entry by entry with mcstastox.

```python
from mcdgs.instruments.trex.model import TRex
from mcdgs.mcstas import read_mcstas_map
from mcdgs.reduction.binning import monitor_bins

bins = monitor_bins(TRex().sample_monitor.bin_limits)  # "toa" in s, "wavelength" in Å
sample_map = read_mcstas_map("sample_monitor", run, bins)  # (toa, wavelength)
spectrum = sample_map.sum("toa")  # the variances add
```

Of a scan, the maps come along the scan variable, with the values of the scan as its coordinate,
in the unit the record of the run keeps for it:

```python
scan_maps = read_mcstas_map("sample_monitor", scan, bins)  # (omega, toa, wavelength)
```

`plot_bender_scan` in [`mcdgs.diagnostics.trex.bender_scan`](../reference/diagnostics.md) reads a
scan this way, and the module is a complete example, from a scan to a figure: see
[Bender](../instruments/trex/bender-scan.md#the-scan).

## Detector pixels

[`read_pixels`](../reference/mcstas.md) reads a detector into [scipp](https://scipp.github.io). The
detector must record an event list with pixel ids, as `banana_detector` of T-Rex does:

```python
Monitor_nD(
    radius=3.0,
    yheight=2.2,
    options=(
        "mantid banana, theta limits=[-36 144] bins=288, y bins=136, "
        "neutron pixel min=0 t lambda, list all neutrons"
    ),
    filename="banana_detector",
    at=[0, -0.2978, 0],
    relative="sample_pos",
)
```

`read_pixels` takes the run folder and the names of the detector, source and sample components, and
returns one `(parameters, pixels)` per entry:

```python
from mcdgs.mcstas import read_pixels

points = read_pixels(run.folder, "banana_detector", "source", "sample_pos")
for parameters, pixels in points:
    print(parameters, pixels.sum().value)
```

`parameters` holds the numeric instrument parameters of the entry. `pixels` holds the event weights
summed per pixel, only for pixels that were hit, with coordinates such as `position`,
`sample_position` and `source_position`.

A run started from a beam file has no source in its file. Give the direction of the beam at the
sample, in the frame of `sample_pos`:

```python
points = read_pixels(
    run.folder, "banana_detector", "source", "sample_pos", beam_direction=[0, 0, 1]
)
```
