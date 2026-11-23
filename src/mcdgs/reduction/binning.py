"""Bin edges: for histograms of reduced data, and those of the map of a monitor.

:func:`bin_edges` gives evenly spaced edges, such as in Q and energy transfer, and
:func:`monitor_bins` the edges of the axes of a McStas monitor, so that anything else is read on
the bins of McStas. :func:`bins_title` writes the bins of a histogram for the title of its
figure.
"""

from collections.abc import Mapping

import numpy as np
import scipp as sc

__all__ = ["MONITOR_COORDINATES", "bin_edges", "bins_title", "monitor_bins"]

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


def bins_title(bins: Mapping[str, sc.Variable], per_line: int = 2) -> str:
    """Write the bins of a histogram as ``name=(start, stop, step) [unit]``, for its title.

    It is the form of the bins themselves, as :func:`bin_edges` takes them: a dimension that is
    binned has its step, one that is summed over, a single bin, has none, ``name=(start, stop)``.
    So the title tells which ranges a figure shows and which it sums over.

    Parameters
    ----------
    bins:
        The bin edges of each dimension, evenly spaced, by the name to write, such as
        ``{"Qx": ..., "Qy": ...}``. Evenly spaced points, such as the angles of a scan, are
        written the same way.
    per_line:
        How many dimensions to write on a line.

    Returns
    -------
    :
        The bins, ``per_line`` to a line, without a unit where there is none.

    Examples
    --------
    >>> print(bins_title({
    ...     "Qx": bin_edges("Qx", -2.0, 1.5, 0.1, unit="1/Å"),
    ...     "Qy": bin_edges("Qy", -0.5, 0.5, unit="1/Å"),
    ...     "E/Ei": bin_edges("E/Ei", -0.8, 0.8),
    ... }))
    Qx=(-2, 1.5, 0.1) [1/Å], Qy=(-0.5, 0.5) [1/Å]
    E/Ei=(-0.8, 0.8)
    """
    parts = []
    for name, edges in bins.items():
        values = edges.values
        numbers = [f"{values[0]:.3g}", f"{values[-1]:.3g}"]
        if values.size > 2:
            numbers.append(f"{values[1] - values[0]:.3g}")
        unit = "" if edges.unit in (None, sc.units.dimensionless) else f" [{edges.unit}]"
        parts.append(f"{name}=({', '.join(numbers)}){unit}")
    return "\n".join(", ".join(parts[i : i + per_line]) for i in range(0, len(parts), per_line))
