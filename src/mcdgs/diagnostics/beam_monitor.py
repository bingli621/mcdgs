"""The count rate of a beam monitor, from a time and wavelength map, for any instrument.

:func:`count_rate` reads from a map of the beam the rate a real monitor at that place would
count, on average and at the peak of the pulse, for an efficiency ε(λ) or for the flux itself.
:func:`print_rates` prints those of the flux, and of an example monitor with ε = C λ. The beam
monitor diagnostic of each instrument gives its own map.
"""

from collections.abc import Callable, Mapping
from dataclasses import dataclass

import numpy as np
import scipp as sc
from numpy.typing import NDArray

from mcdgs.core.components.source import ESS_SOURCE_FREQUENCY

# 1/Å: C of the example efficiency ε = C λ, 1e-4 at 2 Å, within the 1e-5 to 1e-3 of beam
# monitors
EFFICIENCY_PER_ANGSTROM = 5e-5
DEAD_TIME = 100e-9  # s: the dead time of the example monitor


@dataclass(frozen=True)
class CountRate:
    """What a monitor counts, from a time and wavelength map.

    Parameters
    ----------
    average:
        The rate averaged over the pulses, in counts/s.
    peak:
        The highest rate while a pulse passes, in counts/s.
    peak_time:
        When in the pulse it is reached, in ms: the center of its time bin.
    """

    average: float
    peak: float
    peak_time: float


# --8<-- [start:count_rate]
def count_rate(
    intensity: NDArray,
    bins: Mapping[str, sc.Variable],
    efficiency: Callable[[NDArray], NDArray] | None = None,
    frequency: float = ESS_SOURCE_FREQUENCY,
) -> CountRate:
    """Return the average and the peak rate a monitor counts, from a map of the beam on it.

    Each bin of the map is in n/s, averaged over the pulses. Weighted by the efficiency at its
    wavelength, it is in counts/s; summed over wavelength, a time profile. Divided by the
    frequency and by the width of its time bin, a bin of that profile is the rate while the
    pulse passes it.

    Parameters
    ----------
    intensity:
        The map, in n/s per bin, with the shape (time bins, wavelength bins).
    bins:
        Its bin edges, ``"toa"`` and ``"wavelength"``, in that order.
    efficiency:
        The efficiency ε of the monitor at a wavelength in Å, applied to the centers of the
        wavelength bins. If ``None``, 1: the rate of neutrons through the opening.
    frequency:
        The pulse frequency of the source, in Hz.

    Returns
    -------
    :
        The average and the peak rate, and the time of the peak.
    """
    time_edges = bins["toa"].to(unit="s").values
    wavelength_centers = sc.midpoints(bins["wavelength"]).to(unit="Å").values
    counted = intensity if efficiency is None else intensity * efficiency(wavelength_centers)
    profile = counted.sum(axis=1)  # counts/s per time bin, averaged over the pulses
    instantaneous = profile / frequency / np.diff(time_edges)  # counts/s while the pulse passes
    peak_bin = int(np.argmax(instantaneous))
    peak_time = float(sc.midpoints(bins["toa"]).to(unit="ms").values[peak_bin])
    return CountRate(
        average=float(profile.sum()), peak=float(instantaneous[peak_bin]), peak_time=peak_time
    )
    # --8<-- [end:count_rate]


def print_rates(intensity: NDArray, bins: Mapping[str, sc.Variable], area: float) -> None:
    """Print the rates of the flux, and of the example monitor with ε = C λ.

    Parameters
    ----------
    intensity, bins:
        The map and its bin edges, as in :func:`count_rate`.
    area:
        The opening of the monitor, in cm².
    """
    flux = count_rate(intensity, bins)
    print(
        f"flux: {flux.average:.3g} n/s, {flux.average / ESS_SOURCE_FREQUENCY:.3g} n per pulse, "
        f"{flux.average / area:.3g} n/s/cm², peak {flux.peak:.3g} n/s at {flux.peak_time:.2f} ms"
    )
    counted = count_rate(intensity, bins, lambda wavelength: EFFICIENCY_PER_ANGSTROM * wavelength)
    print(
        f"ε = {EFFICIENCY_PER_ANGSTROM:g} λ/Å: {counted.average:.3g} counts/s, "
        f"peak {counted.peak:.3g} counts/s at {counted.peak_time:.2f} ms, "
        f"peak × dead time {counted.peak * DEAD_TIME:.3g}"
    )
