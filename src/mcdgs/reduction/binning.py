"""Bin edges for histograms of reduced data, such as in Q and energy transfer."""

import numpy as np
import scipp as sc

__all__ = ["bin_edges"]


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
