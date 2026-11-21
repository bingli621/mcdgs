"""Read and write McStas NeXus (HDF5) output.

- :func:`read_monitor_map` reads the map of a monitor of one entry, on its bin edges, with h5py.
- :func:`read_scan_maps` reads that map of every point of a scan, along the scan variable.
- :func:`read_events` reads the neutrons of a detector, by pixel, with mcstastox.
- :func:`read_pixels` reads the detector pixels of every scan point, with mcstastox.
- :func:`numeric_parameters` keeps the instrument parameters of a file that are numbers.
- :func:`merge_entries` merges single-run files into one scan file.

Nothing here knows an instrument: the components to read are named by the caller.

A scan file holds one entry per scan point (entry1, entry2, ...). Reading it needs a mcstastox
with ``entry_number`` support, which pyproject.toml pins. ``mcstasscript.load_data`` reads only
a single run, so the maps of a scan are read from here; :func:`read_monitor_map` reads them as
:func:`.read_mcstas_map` does, in n/s per bin on the bin edges of :func:`.monitor_bins`.
"""

from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import h5py
import mcstastox
import numpy as np
import scipp as sc
from numpy.typing import ArrayLike

from mcdgs.mcstas.runs import Run


def numeric_parameters(parameters: Mapping[str, Any]) -> dict[str, float]:
    """Keep the instrument parameters that are numbers.

    Parameters
    ----------
    parameters:
        Instrument parameters by name.

    Returns
    -------
    :
        The numeric ones as floats. Others, e.g. the beam file name, are dropped.
    """
    numeric = {}
    for key, value in parameters.items():
        try:
            numeric[key] = float(value)
        except (TypeError, ValueError):
            pass
    return numeric


def read_monitor_map(
    file: str | Path, name: str, bins: Mapping[str, sc.Variable], entry: int = 1
) -> sc.DataArray:
    """Read the map of a monitor of one entry of a McStas NeXus file, in n/s per bin.

    Parameters
    ----------
    file:
        The McStas NeXus file, such as :attr:`.Run.data_file`.
    name:
        The name of the monitor component. Its group is found by the name after its number in
        beam order, so ``"sample_monitor"`` does not find ``0073_pre_sample_monitor``.
    bins:
        The bin edges of the map, by coordinate, in the order of the axes of the monitor, such as
        ``"toa"`` then ``"wavelength"`` of ``"time ... lambda ..."``: as :func:`.monitor_bins`
        gives them from its ``bin_limits``.
    entry:
        The entry to read, from 1: of a scan file, the scan point.

    Returns
    -------
    :
        The intensity per bin, in n/s, with the square of its error as variances. Its dimensions
        are those of ``bins``, in their order, and its coordinates the bin edges.

    Raises
    ------
    ValueError
        If the entry has not one component of that name, the component has not one map, or the
        map does not have the shape of ``bins``.
    """
    with h5py.File(file) as data:
        components = data[f"entry{entry}/instrument/components"]
        assert isinstance(components, h5py.Group)
        found = [group for key, group in components.items() if key.partition("_")[2] == name]
        if len(found) != 1:
            raise ValueError(f"{file} entry{entry}: {len(found)} components {name!r}, not one")
        output = found[0].get("output", {})
        maps = [group for group in output.values() if "data" in group]
        if len(maps) != 1:
            raise ValueError(f"{file} entry{entry}: {name} has {len(maps)} maps, not one")
        intensity = np.asarray(maps[0]["data"][()], dtype=float)  # axes in the order of bins
        error = np.asarray(maps[0]["errors"][()], dtype=float)
    shape = tuple(len(edges) - 1 for edges in bins.values())
    if intensity.shape != shape:
        raise ValueError(
            f"{name}: McStas has a map of {intensity.shape} bins, the bins are {shape}"
        )
    return sc.DataArray(
        sc.array(dims=list(bins), values=intensity, variances=error**2, unit="1/s"),
        coords=dict(bins),
    )


def read_scan_maps(
    run: Run, name: str, bins: Mapping[str, sc.Variable], unit: str | None = None
) -> sc.DataArray:
    """Read the map of a monitor of every point of a scan, as :func:`read_monitor_map`.

    Parameters
    ----------
    run:
        The scan, as :meth:`.Simulation.scan` returns it or :meth:`.Run.load` reads it.
    name, bins:
        The monitor and its bin edges, as in :func:`read_monitor_map`.
    unit:
        The unit of the scan variable, which the record of the run does not keep.

    Returns
    -------
    :
        The maps along the scan variable, then the dimensions of ``bins``, with the values of the
        variable as its coordinate, in the order of the scan.

    Raises
    ------
    ValueError
        If the run is not a scan, or as :func:`read_monitor_map`.
    """
    if run.variable is None or run.values is None:
        raise ValueError(f"{run.folder}: not a scan")
    maps = [
        read_monitor_map(run.data_file, name, bins, entry)
        for entry in range(1, len(run.values) + 1)
    ]
    scan = sc.concat(maps, run.variable)
    scan.coords[run.variable] = sc.array(dims=[run.variable], values=run.values, unit=unit)
    return scan


def _export_events(
    data: mcstastox.Read,
    folder: str | Path,
    detector: str,
    source: str,
    sample: str,
    beam_direction: ArrayLike | None,
) -> sc.DataArray:
    """Export the events of ``detector``, by pixel, with mcstastox, from an open file.

    Without ``source`` in the file, the source is put 1 m upstream of ``sample`` along
    ``beam_direction``, as :func:`read_pixels` describes.

    Raises
    ------
    ValueError
        If the file has no ``source`` and ``beam_direction`` is ``None``.
    """
    if source in data.get_components():
        return data.export_scipp(source, sample, component_name=detector)["events"]
    # no source in the file, e.g. a run from a beam file
    if beam_direction is None:
        raise ValueError(f"{folder}: no component {source!r}; give beam_direction to read it")
    direction = np.asarray(beam_direction, float)
    direction = direction / np.linalg.norm(direction)
    sample_pos = data.get_global_component_coordinates(sample)
    grouped = data.export_scipp(sample, sample, component_name=detector)["events"]
    grouped.coords["source_position"] = sc.vector(sample_pos - direction, unit="m")  # 1 m upstream
    return grouped


def read_events(
    folder: str | Path,
    detector: str,
    source: str,
    sample: str,
    beam_direction: ArrayLike | None = None,
    entry: int = 1,
) -> sc.DataArray:
    """Read the neutrons of a detector of a McStas NeXus file, by pixel.

    Parameters
    ----------
    folder:
        Run folder with the ``mccode.h5`` to read.
    detector:
        Name of the detector component, which must record an event list.
    source:
        Name of the component taken as the source. Only ``source_position`` comes from it, so
        any component on the incident beam will do: one just before the sample gives the
        direction of the beam there, where a guide bends it.
    sample:
        Name of the component at the sample position.
    beam_direction:
        Direction of the beam at the sample, for a file without ``source``, as in
        :func:`read_pixels`. If ``None``, a file without ``source`` raises.
    entry:
        The entry to read, from 1: of a scan file, the scan point.

    Returns
    -------
    :
        The neutrons, binned by ``pixel_id``, each with its weight in n/s and its time
        ``time_at_detector``, in ms. The coords are those mcstastox gives, in the frame of the file,
        including ``position``, ``sample_position`` and ``source_position``.

    Raises
    ------
    ValueError
        If the file has no ``source`` and ``beam_direction`` is ``None``.
    """
    with mcstastox.Read(str(folder), entry_number=entry) as data:
        binned = _export_events(data, folder, detector, source, sample, beam_direction)
    binned.bins.coords["time_at_detector"] = binned.bins.coords.pop("t").to(unit="ms")
    binned.bins.unit = "1/s"  # McStas weights are n/s
    return binned


def read_pixels(
    folder: str | Path,
    detector: str,
    source: str,
    sample: str,
    beam_direction: ArrayLike | None = None,
) -> list[tuple[dict[str, float], sc.DataArray]]:
    """Read the detector pixels of every scan point of a McStas NeXus file.

    Parameters
    ----------
    folder:
        Run folder with the ``mccode.h5`` to read.
    detector:
        Name of the detector component, which must record an event list.
    source:
        Name of the source component.
    sample:
        Name of the component at the sample position.
    beam_direction:
        Direction of the beam at the sample, for a file without ``source``, such as a run started
        from a beam file at the sample. It is in the frame of the file: for such a run that is
        the frame of ``sample_pos``, not the lab frame (see :class:`.Simulation`). The source is
        then put 1 m upstream of the sample along it, which gives the incident direction; the
        distance does not matter. If ``None``, a file without ``source`` raises.

    Returns
    -------
    :
        One ``(parameters, pixels)`` per scan point:

        - ``parameters``: the numeric instrument parameters of the point.
        - ``pixels``: event weights summed per pixel, over ``pixel_id``. Only pixels that were
          hit are included. The coords are those mcstastox gives, including ``position``,
          ``sample_position`` and ``source_position``.

    Raises
    ------
    ValueError
        If the file has no ``source`` and ``beam_direction`` is ``None``.
    """
    with mcstastox.Read(str(folder)) as data:
        n_entries = data.get_number_of_entries()

    points = []
    for entry in range(1, n_entries + 1):
        with mcstastox.Read(str(folder), entry_number=entry) as data:
            params = numeric_parameters(data.get_instrument_parameters())
            grouped = _export_events(data, folder, detector, source, sample, beam_direction)
        points.append((params, grouped.bins.sum()))  # events are grouped by pixel id
    return points


def merge_entries(files: Iterable[str | Path], merged: Path) -> Path:
    """Merge single-entry McStas files into one scan file.

    Parameters
    ----------
    files:
        McStas NeXus files, each with one entry ``entry1``.
    merged:
        File to write, overwritten if it exists.

    Returns
    -------
    :
        ``merged``, with the entries ``entry1``, ``entry2``, ... in the order of ``files``.
    """
    with h5py.File(merged, "w") as out:
        for n, file in enumerate(files, start=1):
            with h5py.File(file) as part:
                part.copy(part["entry1"], out, name=f"entry{n}")
    return merged
