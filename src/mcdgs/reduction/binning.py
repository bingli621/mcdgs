"""Bin edges: for histograms of reduced data, and those of the map of a monitor.

:func:`bin_edges` gives evenly spaced edges, such as in Q and energy transfer, and
:func:`monitor_bins` the edges of the axes of a McStas monitor, so that anything else is read on
the bins of McStas.
"""

from collections.abc import Mapping

import numpy as np
import scipp as sc

__all__ = ["MONITOR_COORDINATES", "bin_edges", "monitor_bins"]

# the variables of a McStas monitor that are coordinates of the neutrons, by McStas name: the
# coordinate, as scippneutron and tof name it, and the unit of McStas
MONITOR_COORDINATES: dict[str, tuple[str, str]] = {
    "time": ("toa", "s"),
    "lambda": ("wavelength", "Å"),
}


def bin_edges(
    dim: str,
    start: float,
    stop: float,
    step: float | None = None,
    unit: str | None = "dimensionless",
) -> sc.Variable:
    """Return evenly spaced bin edges from ``start`` to ``stop``, both included.

    Parameters
    ----------
    dim:
        The dimension of the edges.
    start, stop:
        The first and the last edge.
    step:
        The width of the bins, which must divide ``stop - start``. If ``None``, a single bin from
        ``start`` to ``stop``, to integrate over the dimension.
    unit:
        The unit of the edges, such as ``"1/Å"`` or ``"meV"``; by default none, as for h, k
        and l.

    Returns
    -------
    :
        The edges, along ``dim``.

    Raises
    ------
    ValueError
        If ``stop`` is not above ``start``, or ``step`` is not positive or does not divide the
        range.

    Examples
    --------
    >>> bin_edges("energy_transfer", -2.0, 2.0, 1.0, unit="meV").values
    array([-2., -1.,  0.,  1.,  2.])
    """
    if not stop > start:
        raise ValueError(f"{dim}: stop must be above start, got {start} to {stop}")
    if step is None:
        bins = 1
    else:
        if not step > 0:
            raise ValueError(f"{dim}: step must be positive, got {step}")
        exact = (stop - start) / step
        bins = round(exact)
        if bins < 1 or not np.isclose(exact, bins, rtol=0.0, atol=1e-9 * max(1.0, exact)):
            raise ValueError(f"{dim}: step {step} does not divide {start} to {stop}")
    return sc.linspace(dim, start, stop, bins + 1, unit=unit)


def monitor_bins(limits: Mapping[str, tuple[float, float, int]]) -> dict[str, sc.Variable]:
    """Return the bin edges of each axis of a McStas monitor, by the coordinate of the neutrons.

    Parameters
    ----------
    limits:
        Start, stop and number of bins of each axis of the monitor, by McStas variable, such as
        ``{"time": (0.01, 0.1, 1800), "lambda": (0.1, 10.0, 495)}``: the ``bin_limits`` of a
        ``Monitor_nD``.

    Returns
    -------
    :
        The edges by coordinate, such as ``"toa"`` in s and ``"wavelength"`` in Å, in the order
        of ``limits``: that of the axes of the map of McStas. An axis that is no coordinate of
        the neutrons, one not in ``MONITOR_COORDINATES``, is left out.

    Examples
    --------
    >>> bins = monitor_bins({"x": (-0.05, 0.05, 10), "lambda": (0.5, 5.0, 50)})
    >>> list(bins), str(bins["wavelength"].unit), bins["wavelength"].sizes
    (['wavelength'], 'Å', {'wavelength': 51})
    """
    bins = {}
    for variable, (start, stop, n_bins) in limits.items():
        if variable in MONITOR_COORDINATES:
            coordinate, unit = MONITOR_COORDINATES[variable]
            bins[coordinate] = sc.linspace(coordinate, start, stop, n_bins + 1, unit=unit)
    return bins
