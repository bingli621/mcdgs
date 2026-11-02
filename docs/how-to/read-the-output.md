# Read the output

Every run writes a McStas NeXus (HDF5) file, `mccode.h5`. A scan has one entry per value, `entry1`,
`entry2` and so on, and a single run has one.

## A monitor

Monitors are read with `h5py`. This reads the `Monitor_nD` with `filename="spectrum"` from
[Install and first run](../getting-started.md):

```python
import h5py

with h5py.File(h5_file) as file:
    components = file["entry1/instrument/components"]
    monitor = next(c for name, c in components.items() if name.endswith("_spectrum"))
    spectrum = monitor["output/spectrum_L"]
    wavelengths = spectrum["Wavelength__Angs_"][()]  # bin centres, in Å
    intensity = spectrum["data"][()]
    error = spectrum["errors"][()]
```

Component groups are numbered in beam order, such as `0067_spectrum`, which is why this finds the
group by the end of its name. McStas adds a suffix for the variable to the file name, `_L` for
wavelength. To see what a monitor wrote, use `list(monitor["output"])`.

To read every entry of a scan, count through the entries, since sorted by name `entry10` comes
before `entry2`:

```python
with h5py.File(scan_file) as file:
    for n in range(1, len(file) + 1):
        entry = file[f"entry{n}"]
```

[`diagnosis/trex/bender_efficiency.py`](https://github.com/bingli621/mcdgs/blob/main/diagnosis/trex/bender_efficiency.py)
is a complete example, from a scan to a figure.

## Detector pixels

[`read_pixels`](../reference/mcstas.md) reads a detector into [scipp](https://scipp.github.io). The
detector must record an event list with pixel ids, as `Small` does in
[Add an instrument](add-an-instrument.md):

```python
Monitor_nD(
    xwidth=0.1,
    yheight=0.1,
    options="x bins=4, y bins=3, list all, pixel id, t",
    at=[0, 0, 1],
    relative="sample_pos",
)
```

`read_pixels` takes the run folder and the names of the detector, source and sample components, and
returns one `(parameters, pixels)` per entry:

```python
from mcdgs.mcstas.nexus_io import read_pixels

points = read_pixels(h5_file.parent, "detector", "source", "sample_pos")
for parameters, pixels in points:
    print(parameters, pixels.sum().value)
```

`parameters` holds the numeric instrument parameters of the entry. `pixels` holds the event weights
summed per pixel, only for pixels that were hit, with coordinates such as `position`,
`sample_position` and `source_position`.

A run started from a beam file has no source in its file. Give the direction of the beam at the
sample, in the frame of `sample_pos`:

```python
points = read_pixels(h5_file.parent, "detector", "source", "sample_pos", beam_direction=[0, 0, 1])
```
