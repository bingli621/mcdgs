"""The solid angle of rectangular detector pixels, seen from the sample.

:func:`rectangle_solid_angle` takes the center and size of each pixel, and the way the detector
faces the sample: a cylinder around a vertical axis through the sample, as the banana detector of
T-Rex, or a sphere around it. The solid angle follows from the four corners of the pixel, split
into two triangles, each with the formula of A. Van Oosterom and J. Strackee, "The Solid Angle of
a Plane Triangle", IEEE Trans. Biomed. Eng. 30, 125 (1983).
"""

from typing import Literal

import numpy as np
import scipp as sc
from numpy.typing import NDArray

__all__ = ["DetectorShape", "pixel_normal", "rectangle_solid_angle"]

DetectorShape = Literal["cylinder", "sphere"]


def pixel_normal(relative_position: sc.Variable, shape: DetectorShape) -> sc.Variable:
    """Return the outward normal of each pixel, not normalized.

    Parameters
    ----------
    relative_position:
        The center of each pixel, from the sample.
    shape:
        ``"cylinder"``, around a vertical (y) axis through the sample: the normal is horizontal;
        or ``"sphere"``, around the sample: the normal points away from it.

    Returns
    -------
    :
        The normal of each pixel, in the unit of ``relative_position``.

    Raises
    ------
    ValueError
        If ``shape`` is neither.
    """
    match shape:
        case "cylinder":
            return sc.spatial.as_vectors(
                relative_position.fields.x,
                sc.zeros_like(relative_position.fields.y),
                relative_position.fields.z,
            )
        case "sphere":
            return relative_position.copy()
        case _:
            raise ValueError(f"Unknown detector shape {shape!r}: 'cylinder' or 'sphere'")


def _unit(vectors: sc.Variable) -> sc.Variable:
    """Return the vectors divided by their length, without unit."""
    return vectors / sc.norm(vectors)


def _triangle_solid_angle(a: sc.Variable, b: sc.Variable, c: sc.Variable) -> NDArray[np.float64]:
    """Return the solid angle of the triangle of the unit vectors a, b and c.

    It is the formula of Van Oosterom and Strackee, ``tan(Ω / 2) = |a · (b × c)| / (1 + a · b
    + b · c + c · a)``, which keeps its precision for small triangles.
    """
    numerator = np.abs(sc.dot(a, sc.cross(b, c)).values)
    denominator = 1 + (sc.dot(a, b) + sc.dot(b, c) + sc.dot(c, a)).values
    return 2 * np.arctan2(numerator, denominator)


def rectangle_solid_angle(
    position: sc.Variable,
    sample_position: sc.Variable,
    width: sc.Variable,
    height: sc.Variable,
    shape: DetectorShape = "cylinder",
) -> sc.Variable:
    """Return the solid angle of rectangular pixels, seen from the sample.

    Each pixel is a flat rectangle, centered on its position and facing along
    :func:`pixel_normal`: its width is horizontal and its height across it. The solid angle is
    exact for such a rectangle, as long as it is below 2π.

    Parameters
    ----------
    position:
        The center of each pixel, as vectors.
    sample_position:
        Where the sample is.
    width, height:
        The size of each pixel. They may be scalars, or have the dimensions of ``position``.
    shape:
        How the detector faces the sample, as :func:`pixel_normal`.

    Returns
    -------
    :
        The solid angle of each pixel, in sr, with the dimensions of ``position``.

    Examples
    --------
    A square of 2 m at 1 m, facing the sample, spans 4 arcsin(1/2) = 2π/3 sr:

    >>> omega = rectangle_solid_angle(
    ...     sc.vector([0.0, 0.0, 1.0], unit="m"),
    ...     sc.vector([0.0, 0.0, 0.0], unit="m"),
    ...     sc.scalar(2.0, unit="m"),
    ...     sc.scalar(2.0, unit="m"),
    ... )
    >>> print(f"{omega.value:.5f} {omega.unit}")
    2.09440 sr
    """
    relative = position - sample_position.to(unit=position.unit)
    normal = pixel_normal(relative, shape)
    across = _unit(  # horizontal, across the normal
        sc.spatial.as_vectors(normal.fields.z, sc.zeros_like(normal.fields.y), -normal.fields.x)
    )
    up = _unit(sc.cross(normal, across))
    half_width = across * (width.to(unit=position.unit) / 2)
    half_height = up * (height.to(unit=position.unit) / 2)

    # the corners, around the rectangle, as unit vectors from the sample
    c1 = _unit(relative - half_width + half_height)
    c2 = _unit(relative + half_width + half_height)
    c3 = _unit(relative + half_width - half_height)
    c4 = _unit(relative - half_width - half_height)
    # two triangles, split along the diagonal from c1 to c3
    omega = _triangle_solid_angle(c1, c2, c3) + _triangle_solid_angle(c1, c3, c4)
    return sc.array(dims=relative.dims, values=omega, unit="sr")
