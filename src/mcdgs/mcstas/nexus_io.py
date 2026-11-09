"""Read and write McStas NeXus (HDF5) output.

- :func:`read_events` reads the neutrons of a detector, by pixel, with mcstastox.
- :func:`read_pixels` reads the detector pixels of every scan point, with mcstastox.
- :func:`merge_entries` merges single-run files into one scan file.

Nothing here knows an instrument: the components to read are named by the caller.

A scan file holds one entry per scan point (entry1, entry2, ...). Reading it needs a mcstastox
with ``entry_number`` support, which pixi.toml and pyproject.toml pin.
"""

from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import h5py
import mcstastox
import numpy as np
import scipp as sc
from numpy.typing import ArrayLike


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
    try:
        return data.export_scipp(source, sample, component_name=detector)["events"]
    except ValueError as error:  # no source in the file, e.g. a run from a beam file
        if beam_direction is None:
            raise ValueError(
                f"{folder}: no component {source!r}; give beam_direction to read it"
            ) from error
        sample_pos = data.get_global_component_coordinates(sample)
        upstream = np.asarray(beam_direction, float)
        upstream = upstream / np.linalg.norm(upstream)
        grouped = data.export_scipp(sample, sample, component_name=detector)["events"]
        grouped.coords["source_position"] = sc.vector(sample_pos - upstream, unit="m")
        return grouped


def read_events(
    folder: str | Path,
    detector: str,
    source: str,
    sample: str,
    beam_direction: ArrayLike | None = None,
) -> sc.DataArray:
    """Read the neutrons of a detector of a McStas NeXus file, by pixel.

    Parameters
    ----------
    folder:
        Run folder with the ``mccode.h5`` to read. Of a scan file, the first scan point.
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
    with mcstastox.Read(str(folder)) as data:
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
    points = []
    with mcstastox.Read(str(folder)) as data:
        n_entries = data.get_number_of_entries()

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
