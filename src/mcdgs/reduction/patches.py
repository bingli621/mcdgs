"""The patches of a time and wavelength map, and the incident energy of the sub-pulses among them.

A monitor of a repetition-rate multiplication instrument sees the sub-pulses of each pulse as
patches of its time and wavelength map, a few bins wide and some ms apart, and, much weaker,
the patches of the neutrons that leak between them. :func:`find_patches` gives the center of
mass, width and intensity of each patch. :func:`subpulse_flight` keeps the patches strong enough
to be sub-pulses at two monitors, pairs them, and gives each sub-pulse the incident energy of its
flight between the monitors, with :func:`.incident_flight`.
"""

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np
import scipp as sc
from numpy.typing import NDArray

from mcdgs.reduction.incident_energy import incident_flight

__all__ = ["MIN_FRACTION", "PATCH_GAP", "Patches", "find_patches", "subpulse_flight"]

PATCH_GAP = 1.0  # ms: a longer gap in time starts a new patch; they are at least 2 ms apart
MIN_FRACTION = 0.02  # a patch weaker than this fraction of the strongest is not a sub-pulse


@dataclass(frozen=True)
class Patches:
    """The center of mass of each patch of a time and wavelength map, and its width.

    Each field has one value per patch, in time order. A width is the RMS of the patch about its
    center, with the spread within a bin, its width over √12, added in quadrature: a patch one
    bin wide is that wide, not a point.

    Parameters
    ----------
    time, time_width:
        The mean time of arrival and its width, in ms.
    wavelength, wavelength_width:
        The mean wavelength and its width, in Å.
    intensity:
        The intensity of the patch, in n/s.
    """

    time: NDArray
    time_width: NDArray
    wavelength: NDArray
    wavelength_width: NDArray
    intensity: NDArray

    def strongest(self, min_fraction: float = MIN_FRACTION) -> "Patches":
        """Return the patches with at least ``min_fraction`` of the intensity of the strongest.

        Examples
        --------
        >>> patches = Patches(*np.array([[1.0, 2.0, 3.0]] * 4), np.array([100.0, 1.0, 50.0]))
        >>> patches.strongest(0.02).intensity.tolist()
        [100.0, 50.0]
        """
        keep = self.intensity >= min_fraction * self.intensity.max()
        return Patches(
            self.time[keep],
            self.time_width[keep],
            self.wavelength[keep],
            self.wavelength_width[keep],
            self.intensity[keep],
        )


def _mean_and_width(weights: NDArray, edges: NDArray) -> tuple[float, float]:
    """Return the weighted mean of the bin centers and the width, with the spread within a bin."""
    centers = (edges[:-1] + edges[1:]) / 2
    mean = float(np.average(centers, weights=weights))
    variance = np.average((centers - mean) ** 2, weights=weights) + np.diff(edges)[0] ** 2 / 12
    return mean, float(np.sqrt(variance))


def find_patches(
    intensity: NDArray, bins: Mapping[str, sc.Variable], gap: float = PATCH_GAP
) -> Patches:
    """Find the patches of a time and wavelength map: the runs of time bins with neutrons.

    Parameters
    ----------
    intensity:
        The map, with the shape (time bins, wavelength bins).
    bins:
        Its bin edges, ``"toa"`` and ``"wavelength"``, evenly spaced.
    gap:
        The longest gap in time within a patch, in ms; a longer one starts a new patch.

    Returns
    -------
    :
        The center of mass of each patch, its width and its intensity.

    Examples
    --------
    >>> bins = {
    ...     "toa": sc.linspace("toa", 0.0, 6.0, 7, unit="ms"),
    ...     "wavelength": sc.linspace("wavelength", 1.0, 3.0, 3, unit="Å"),
    ... }
    >>> intensity = np.array([[1, 0], [1, 0], [0, 0], [0, 0], [0, 2], [0, 0]])
    >>> patches = find_patches(intensity, bins, gap=1.0)
    >>> patches.time.tolist(), patches.wavelength.tolist(), patches.intensity.tolist()
    ([1.0, 4.5], [1.5, 2.5], [2.0, 2.0])
    """
    times = bins["toa"].to(unit="ms").values
    wavelengths = bins["wavelength"].to(unit="Å").values
    filled = np.flatnonzero(intensity.sum(axis=1) > 0)
    new_patch = np.flatnonzero(np.diff(times[filled]) > gap + np.diff(times)[0]) + 1
    found = []
    for run in np.split(filled, new_patch) if filled.size else []:
        patch = intensity[run[0] : run[-1] + 1]
        time = _mean_and_width(patch.sum(axis=1), times[run[0] : run[-1] + 2])
        wavelength = _mean_and_width(patch.sum(axis=0), wavelengths)
        found.append((*time, *wavelength, float(patch.sum())))
    columns = np.array(found, dtype=float).reshape(-1, 5).T
    return Patches(*columns)


def subpulse_flight(
    start: Patches,
    start_position: sc.Variable,
    end: Patches,
    end_position: sc.Variable,
    min_fraction: float = MIN_FRACTION,
) -> sc.DataArray:
    """Pair the sub-pulses at two monitors, and compute their flight between them.

    The sub-pulses are the patches with at least ``min_fraction`` of the intensity of the
    strongest, as :meth:`Patches.strongest`. The flight between the monitors must be shorter
    than the time between sub-pulses, so that the n-th sub-pulse at one is the n-th at the
    other.

    Parameters
    ----------
    start, end:
        The patches at the monitor where the flight starts, and where it ends.
    start_position, end_position:
        Where the monitors are, as vectors.
    min_fraction:
        The weakest patch that counts as a sub-pulse, as a fraction of the strongest.

    Returns
    -------
    :
        The intensity of each sub-pulse at the start, in n/s, along ``subpulse``, with
        ``time_at_start`` and ``time_at_end``, in ms, ``wavelength_at_start``, its mean
        wavelength there, and the coordinates of :func:`.incident_flight`: ``Ltotal``, ``tof``,
        ``incident_wavelength`` and ``incident_energy``.

    Raises
    ------
    ValueError
        If the monitors do not have as many sub-pulses, or one reaches the end first.

    Examples
    --------
    A sub-pulse of 2 Å takes 1.5167 ms over 3 m; the weak patch between is left out:

    >>> def patches(time: list[float], intensity: list[float]) -> Patches:
    ...     ones = np.ones(len(time))
    ...     return Patches(np.array(time), ones, 2 * ones, ones, np.array(intensity))
    >>> flight = subpulse_flight(
    ...     patches([10.0, 12.0], [100.0, 1.0]),
    ...     sc.vector([0.0, 0.0, 0.0], unit="m"),
    ...     patches([11.5167], [90.0]),
    ...     sc.vector([0.0, 0.0, 3.0], unit="m"),
    ... )
    >>> print(f"{flight.coords['incident_wavelength'].to(unit='Å').values[0]:.3f}")
    2.000
    """
    start, end = start.strongest(min_fraction), end.strongest(min_fraction)
    if start.time.size != end.time.size:
        raise ValueError(f"{start.time.size} sub-pulses at the start, {end.time.size} at the end")
    time_at_start = sc.array(dims=["subpulse"], values=start.time, unit="ms")
    time_at_end = sc.array(dims=["subpulse"], values=end.time, unit="ms")
    flight = incident_flight(start_position, time_at_start, end_position, time_at_end)
    if (flight["tof"] <= sc.scalar(0.0, unit="ms")).any().value:
        raise ValueError("A sub-pulse reaches the end first")
    return sc.DataArray(
        sc.array(dims=["subpulse"], values=start.intensity, unit="1/s"),  # n/s
        coords={
            "time_at_start": time_at_start,
            "time_at_end": time_at_end,
            "wavelength_at_start": sc.array(dims=["subpulse"], values=start.wavelength, unit="Å"),
            **flight,
        },
    )
