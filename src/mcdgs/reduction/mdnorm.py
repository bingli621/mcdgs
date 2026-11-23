"""The coverage of each bin of (|Q|, E): what normalizes a histogram to S(Q, E).

A pixel at the scattering angle 2θ sees, for a sub-pulse of incident energy Ei, every energy
transfer E of the arrival window, each at its own |Q|, from :func:`momentum_transfer`: a line
through (|Q|, E), its trajectory. :func:`q_e_coverage` sums, for each bin, the solid angle of
every pixel times the length in E of its trajectory inside the bin. A histogram of the fraction
of the incident neutrons scattered into each bin, divided by that coverage, is S(Q, E) per
steradian and per meV, the scattering of the sample whatever the part of the detector that sees
each bin. It is the normalization of MDNorm, for |Q| rather than the vector Q. A bin no
trajectory crosses has no coverage: the sub-pulse does not reach it.
"""

import numpy as np
import scipp as sc
from numpy.typing import ArrayLike, NDArray
from scipp import constants

__all__ = [
    "HBAR2_OVER_2M",
    "SAMPLES_PER_BIN",
    "TWO_THETA_STEP",
    "momentum_transfer",
    "q_e_coverage",
]

# meV Å², ħ²/2m of the neutron: E = HBAR2_OVER_2M k²
HBAR2_OVER_2M = float((constants.hbar**2 / (2 * constants.neutron_mass)).to(unit="meV*Å**2").value)
SAMPLES_PER_BIN = 4  # a trajectory is followed at this many points in each bin of E
TWO_THETA_STEP = 0.02  # deg: pixels closer in 2θ than this share one trajectory
E_BLOCK = 256  # the points in E followed at once, to bound the memory


def momentum_transfer(
    incident_energy: float, energy_transfer: ArrayLike, two_theta: ArrayLike
) -> NDArray[np.float64]:
    """Return |Q| of a neutron of an incident energy and energy transfer, scattered by 2θ.

    ``Q² = ki² + kf² - 2 ki kf cos 2θ``, with ``E = HBAR2_OVER_2M k²``. The energy transfer and
    the angle broadcast against each other.

    Parameters
    ----------
    incident_energy:
        Ei, in meV.
    energy_transfer:
        E = Ei - Ef, in meV, below Ei.
    two_theta:
        The scattering angle, in deg.

    Returns
    -------
    :
        |Q|, in 1/Å.

    Examples
    --------
    Elastic and scattered straight back, Q is 2 ki:

    >>> ki = np.sqrt(20.0 / HBAR2_OVER_2M)
    >>> bool(np.isclose(momentum_transfer(20.0, 0.0, 180.0), 2 * ki))
    True
    """
    ki = np.sqrt(incident_energy / HBAR2_OVER_2M)
    kf = np.sqrt((incident_energy - np.asarray(energy_transfer, dtype=float)) / HBAR2_OVER_2M)
    cosine = np.cos(np.radians(np.asarray(two_theta, dtype=float)))
    return np.asarray(np.sqrt(ki**2 + kf**2 - 2 * ki * kf * cosine), dtype=np.float64)


def _trajectories(
    two_theta: NDArray[np.float64], solid_angle: NDArray[np.float64], step: float
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Group the pixels by 2θ, ``step`` wide: the mean 2θ and the summed solid angle of each."""
    group = np.floor(two_theta / step).astype(np.int64)
    groups, index = np.unique(group, return_inverse=True)
    total = np.bincount(index, solid_angle, minlength=groups.size)
    mean = np.bincount(index, solid_angle * two_theta, minlength=groups.size) / total
    return mean, total


def q_e_coverage(
    two_theta: sc.Variable,
    solid_angle: sc.Variable,
    incident_energy: sc.Variable,
    energy_edges: sc.Variable,
    q_edges: sc.Variable,
    window: tuple[float, float],
    samples_per_bin: int = SAMPLES_PER_BIN,
    two_theta_step: float = TWO_THETA_STEP,
) -> sc.Variable:
    """Return the coverage of each bin of (E, |Q|) for one sub-pulse: Σ solid angle × ΔE.

    Each bin of E, cut to the window, ``window`` times Ei, is followed at ``samples_per_bin``
    evenly spaced points; at each, every trajectory adds its solid angle times the spacing of
    the points to the bin of |Q| it reaches there. The lengths in E are exact, the bins of |Q|
    as fine as the points. The pixels are grouped by 2θ, ``two_theta_step``
    wide, into one trajectory each.

    Parameters
    ----------
    two_theta, solid_angle:
        The scattering angle and the solid angle of each pixel, along one dimension.
    incident_energy:
        Ei of the sub-pulse, a scalar.
    energy_edges, q_edges:
        The bin edges in E and in |Q|, evenly spaced.
    window:
        The smallest and the largest ``E / Ei`` of the arrival window.
    samples_per_bin:
        The points each bin of E is followed at.
    two_theta_step:
        The width in 2θ, in deg, of a group of pixels.

    Returns
    -------
    :
        The coverage of each bin, along the dimensions of ``energy_edges`` and ``q_edges``, in
        the unit of ``solid_angle`` times that of ``energy_edges``.

    Raises
    ------
    ValueError
        If the bins in |Q| are not evenly spaced.

    Examples
    --------
    One pixel of 1 sr at 90°, an elastic window: it covers the bin of E = 0, 1 meV wide, over
    the length of the window inside it, 0.2 meV, at the |Q| of 90°:

    >>> coverage = q_e_coverage(
    ...     sc.array(dims=["pixel"], values=[90.0], unit="deg"),
    ...     sc.array(dims=["pixel"], values=[1.0], unit="sr"),
    ...     sc.scalar(20.0, unit="meV"),
    ...     sc.array(dims=["energy_transfer"], values=[-0.5, 0.5], unit="meV"),
    ...     sc.linspace("Q", 0.0, 5.0, 51, unit="1/Å"),
    ...     window=(-0.005, 0.005),
    ... )
    >>> round(float(coverage.sum().value), 6), str(coverage.unit)
    (0.2, 'meV*sr')
    """
    energy = energy_edges.values
    q = q_edges.to(unit="1/Å").values
    q_step = np.diff(q)
    if not np.allclose(q_step, q_step[0]):
        raise ValueError("the bins in |Q| must be evenly spaced")
    ei = float(incident_energy.to(unit=energy_edges.unit).value)
    angles, weights = _trajectories(
        two_theta.to(unit="deg").values, solid_angle.values, two_theta_step
    )
    # each bin of E, cut to the window, is followed at evenly spaced points: the bin each point
    # is in, and the length of E each stands for, so that the lengths add up exactly
    low = np.maximum(energy[:-1], window[0] * ei)
    high = np.minimum(energy[1:], window[1] * ei)
    (inside,) = np.nonzero(high > low)
    fraction = (np.arange(samples_per_bin) + 0.5) / samples_per_bin
    width = high[inside] - low[inside]
    points = (low[inside, None] + width[:, None] * fraction).ravel()
    rows = np.repeat(inside, samples_per_bin)
    spacing = np.repeat(width / samples_per_bin, samples_per_bin)

    n_q = q.size - 1
    total = np.zeros((energy.size - 1) * n_q)
    for start in range(0, points.size, E_BLOCK):
        block = slice(start, start + E_BLOCK)
        reached = momentum_transfer(ei, points[block, None], angles[None, :])
        column = np.floor((reached - q[0]) / q_step[0]).astype(np.int64)
        valid = (column >= 0) & (column < n_q)
        cell = (rows[block, None] * n_q + column)[valid]
        weight = (spacing[block, None] * weights[None, :])[valid]
        total += np.bincount(cell, weight, minlength=total.size)
    values = sc.array(dims=[*energy_edges.dims, *q_edges.dims], values=total.reshape(-1, n_q))
    # in the unit of the solid angle times that of E
    return values * sc.scalar(1.0, unit=solid_angle.unit) * sc.scalar(1.0, unit=energy_edges.unit)
