"""The rays of a beam in a window of position and divergence, and their brilliance.

The rays are those that :meth:`.McStasScriptDiagnostics.rays` reads at a point: the
``RAY_VARIABLES`` of every ray, in the frame of the point. A :class:`Window` is an area across
the beam and a divergence of ± its half divergence from the axis, both horizontally and
vertically.

- :func:`window_intensity` sums the weight of the rays in a window, in each wavelength bin;
- :func:`brilliance` divides it by the area, the solid angle and the width of the bin, in
  n/s/cm²/sr/Å;
- :func:`brilliance_transfer` divides a brilliance by a reference, such as that of the
  moderator, in each bin, and :func:`band_transfer` over a band of wavelengths.
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

# the variables of a ray: position in m, divergence in deg, wavelength in Å and weight in n/s
RAY_VARIABLES = ("x", "y", "dx", "dy", "l", "p")

Rays = Mapping[str, NDArray[np.float64]]  # the RAY_VARIABLES of every ray at one point


@dataclass(frozen=True)
class Window:
    """An area across the beam, and a divergence from its axis, both ways.

    Parameters
    ----------
    half_width, half_height:
        Half the size of the area, in m.
    half_divergence:
        The largest angle from the axis, horizontally and vertically, in deg.
    x_center:
        x of the center of the area, in m.
    """

    half_width: float
    half_height: float
    half_divergence: float
    x_center: float = 0.0

    @property
    def area(self) -> float:
        """The area, in cm²."""
        return 4 * self.half_width * self.half_height * 1e4

    @property
    def solid_angle(self) -> float:
        """The solid angle of the divergence, in sr: a square pyramid, 4 asin(sin² a).

        Examples
        --------
        >>> round(Window(0.005, 0.005, 0.5).solid_angle, 9)
        0.00030461
        """
        return 4 * math.asin(math.sin(math.radians(self.half_divergence)) ** 2)

    def contains(self, rays: Rays) -> NDArray[np.bool_]:
        """Tell which rays are in the window, in its area and its divergence."""
        return (
            (np.abs(rays["x"] - self.x_center) <= self.half_width)
            & (np.abs(rays["y"]) <= self.half_height)
            & (np.abs(rays["dx"]) <= self.half_divergence)
            & (np.abs(rays["dy"]) <= self.half_divergence)
        )


def window_intensity(
    rays: Rays, window: Window, bins: NDArray[np.float64]
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Sum the weight of the rays in a window, in each wavelength bin.

    Parameters
    ----------
    rays:
        The rays at the point, with the ``RAY_VARIABLES``, in the frame of the point.
    window:
        The window.
    bins:
        The edges of the wavelength bins, in Å.

    Returns
    -------
    :
        The intensity, in n/s, and the number of rays in the window, in each bin.

    Examples
    --------
    >>> rays = {"x": np.zeros(3), "y": np.zeros(3), "dx": np.array([0.0, 0.1, 2.0]),
    ...         "dy": np.zeros(3), "l": np.array([1.5, 1.5, 2.5]), "p": np.array([1.0, 2.0, 4.0])}
    >>> intensity, n = window_intensity(rays, Window(0.005, 0.005, 0.5), np.array([1.0, 2.0, 3.0]))
    >>> intensity.tolist(), n.tolist()
    ([3.0, 0.0], [2.0, 0.0])
    """
    inside = window.contains(rays)
    intensity, _ = np.histogram(rays["l"][inside], bins=bins, weights=rays["p"][inside])
    counts, _ = np.histogram(rays["l"][inside], bins=bins)
    return intensity.astype(float), counts.astype(float)


def brilliance(
    rays: Rays, window: Window, bins: NDArray[np.float64]
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Compute the brilliance in a window, in each wavelength bin.

    Parameters
    ----------
    rays:
        The rays at the point, with the ``RAY_VARIABLES``, in the frame of the point.
    window:
        The window.
    bins:
        The edges of the wavelength bins, in Å.

    Returns
    -------
    :
        The brilliance, in n/s/cm²/sr/Å, and the number of rays in the window, in each bin.

    Examples
    --------
    >>> rays = {"x": np.zeros(3), "y": np.zeros(3), "dx": np.array([0.0, 0.1, 2.0]),
    ...         "dy": np.zeros(3), "l": np.array([1.5, 1.5, 1.5]), "p": np.ones(3)}
    >>> window = Window(0.005, 0.005, 0.5)
    >>> b, n = brilliance(rays, window, np.array([1.0, 2.0]))
    >>> n.tolist(), bool(np.isclose(b[0], 2 / window.area / window.solid_angle))
    ([2.0], True)
    """
    intensity, counts = window_intensity(rays, window, bins)
    return intensity / (window.area * window.solid_angle * np.diff(bins)), counts


def brilliance_transfer(
    value: NDArray[np.float64],
    rays: NDArray[np.float64],
    reference: NDArray[np.float64],
    reference_rays: NDArray[np.float64],
    min_rays: int = 20,
) -> NDArray[np.float64]:
    """Divide a brilliance by a reference brilliance, such as that of the moderator, in each bin.

    Parameters
    ----------
    value, rays:
        The brilliance in a window, and its number of rays, from :func:`brilliance`.
    reference, reference_rays:
        The reference brilliance, and its number of rays, on the same bins.
    min_rays:
        The fewest rays a bin needs, in both, to be read.

    Returns
    -------
    :
        The ratio; NaN where either has fewer than ``min_rays`` rays, or the reference is 0.

    Examples
    --------
    >>> brilliance_transfer(np.array([1.0, 1.0]), np.array([50.0, 5.0]),
    ...                     np.array([4.0, 4.0]), np.array([100.0, 100.0]))
    array([0.25,  nan])
    """
    enough = (rays >= min_rays) & (reference_rays >= min_rays) & (reference > 0)
    return np.divide(value, reference, out=np.full_like(value, np.nan), where=enough)


def band_transfer(
    value: NDArray[np.float64],
    reference: NDArray[np.float64],
    bins: NDArray[np.float64],
    band: tuple[float, float],
) -> float:
    """Divide the brilliance in a band of wavelengths by a reference brilliance.

    Parameters
    ----------
    value, reference:
        The brilliance in a window and the reference brilliance, in each bin.
    bins:
        The edges of the bins, in Å.
    band:
        The band, in Å: the bins whose centers are in it.

    Returns
    -------
    :
        The ratio of their sums over the band; NaN if the reference has nothing there.
    """
    centers = (bins[:-1] + bins[1:]) / 2
    in_band = (centers > band[0]) & (centers < band[1])
    total = float(reference[in_band].sum())
    return float(value[in_band].sum()) / total if total > 0 else math.nan
