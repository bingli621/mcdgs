"""What the chopper cascade of T-Rex passes to a detector, from a tof run.

The P- and M-Chopper cut the band of the bandwidth choppers, so the neutrons arrive in
sub-pulses. :func:`count_groups` tells them apart by the gaps between arrival times, and
:func:`summarize` gives how many neutrons reach a detector, which wavelengths and when, as an
:class:`ArrivalSummary`.
"""

from dataclasses import dataclass

import numpy as np
import tof
from numpy.typing import ArrayLike


@dataclass(frozen=True)
class ArrivalSummary:
    """What reaches a detector in one pulse of a tof run.

    Parameters
    ----------
    neutrons:
        How many neutrons reach it.
    wavelength_min, wavelength_max:
        Their shortest and longest wavelength, in Å.
    sub_pulses:
        How many groups they arrive in.
    first_arrival, last_arrival:
        When the first and the last arrive, in ms after the start of the pulse.
    """

    neutrons: int
    wavelength_min: float
    wavelength_max: float
    sub_pulses: int
    first_arrival: float
    last_arrival: float


def count_groups(times: ArrayLike, gap: float) -> int:
    """Count the groups of times, a new one starting after each gap longer than ``gap``.

    Parameters
    ----------
    times:
        The times, in any order. There must be at least one.
    gap:
        The longest gap within a group, in the unit of ``times``.

    Returns
    -------
    :
        The number of groups.

    Examples
    --------
    >>> count_groups([0.0, 0.1, 5.0, 5.2, 10.0], gap=1.0)
    3
    """
    ordered = np.sort(np.asarray(times, dtype=float))
    if ordered.size == 0:
        raise ValueError("there are no times to group")
    return int(np.count_nonzero(np.diff(ordered) > gap)) + 1


def summarize(result: tof.Result, detector: str, gap: float = 1.0) -> ArrivalSummary:
    """Summarize the neutrons of the first pulse that reach a detector.

    Parameters
    ----------
    result:
        The result of a tof run.
    detector:
        The name of the detector.
    gap:
        The longest gap between arrivals within a sub-pulse, in ms. The default is chosen for
        T-Rex: its sub-pulses are about 4 ms apart, one turn of the M-Chopper, while within
        one the gaps are under 0.01 ms.

    Returns
    -------
    :
        Their number, wavelength band, sub-pulses and arrival times.
    """
    neutrons = result[detector].data["pulse", 0]
    reached = neutrons[~neutrons.masks["blocked_by_others"]]
    if reached.size == 0:
        raise ValueError(f"no neutron reaches {detector!r}")
    wavelengths = reached.coords["wavelength"].to(unit="angstrom").values
    arrivals = reached.coords["toa"].to(unit="ms").values
    return ArrivalSummary(
        neutrons=int(reached.size),
        wavelength_min=float(wavelengths.min()),
        wavelength_max=float(wavelengths.max()),
        sub_pulses=count_groups(arrivals, gap),
        first_arrival=float(arrivals.min()),
        last_arrival=float(arrivals.max()),
    )
