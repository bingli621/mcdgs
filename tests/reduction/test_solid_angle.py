import math

import pytest
import scipp as sc

from mcdgs.reduction.solid_angle import pixel_normal, rectangle_solid_angle

ORIGIN = sc.vector([0.0, 0.0, 0.0], unit="m")


def _centered_square(side: float, distance: float) -> float:
    """The solid angle of a square of ``side`` at ``distance``, facing it on its center."""
    a = side / 2
    return 4 * math.asin(a**2 / (a**2 + distance**2))


@pytest.mark.parametrize("side, distance", [(2.0, 1.0), (0.1, 3.0), (0.01, 3.0)])
def test_a_centered_square_has_its_exact_solid_angle(side, distance):
    omega = rectangle_solid_angle(
        sc.vector([0.0, 0.0, distance], unit="m"),
        ORIGIN,
        sc.scalar(side, unit="m"),
        sc.scalar(side, unit="m"),
    )
    assert omega.unit == "sr"
    assert omega.value == pytest.approx(_centered_square(side, distance), rel=1e-9)


def test_a_small_pixel_spans_its_area_over_the_distance_squared():
    # a pixel of 25 mm by 16 mm, 3 m away at 60° in the horizontal plane
    position = sc.vector([3.0 * math.sin(math.pi / 3), 0.0, 3.0 * math.cos(math.pi / 3)], unit="m")
    omega = rectangle_solid_angle(
        position, ORIGIN, sc.scalar(25.0, unit="mm"), sc.scalar(16.0, unit="mm")
    )
    assert omega.value == pytest.approx(0.025 * 0.016 / 9.0, rel=1e-4)


def test_a_cylinder_pixel_above_the_beam_is_seen_at_a_slant():
    # facing horizontally, a pixel 1 m above the beam at 3 m is seen at cos(θ) / r²
    position = sc.vector([0.0, 1.0, 3.0], unit="m")
    size = sc.scalar(0.01, unit="m")
    r2 = 10.0
    cylinder = rectangle_solid_angle(position, ORIGIN, size, size, shape="cylinder")
    sphere = rectangle_solid_angle(position, ORIGIN, size, size, shape="sphere")
    assert cylinder.value == pytest.approx(1e-4 * (3.0 / math.sqrt(r2)) / r2, rel=1e-4)
    assert sphere.value == pytest.approx(1e-4 / r2, rel=1e-4)


def test_the_sample_need_not_be_at_the_origin():
    sample = sc.vector([1.0, 2.0, 160.0], unit="m")
    omega = rectangle_solid_angle(
        sample + sc.vector([0.0, 0.0, 1.0], unit="m"),
        sample,
        sc.scalar(2.0, unit="m"),
        sc.scalar(2.0, unit="m"),
    )
    assert omega.value == pytest.approx(_centered_square(2.0, 1.0), rel=1e-9)


def test_one_solid_angle_per_pixel():
    positions = sc.vectors(dims=["pixel"], values=[[0.0, 0.0, 3.0], [3.0, 0.0, 0.0]], unit="m")
    omega = rectangle_solid_angle(
        positions, ORIGIN, sc.scalar(0.1, unit="m"), sc.scalar(0.1, unit="m")
    )
    assert omega.dims == ("pixel",)
    assert omega.values == pytest.approx([_centered_square(0.1, 3.0)] * 2, rel=1e-9)


def test_an_unknown_shape_raises():
    with pytest.raises(ValueError, match="Unknown detector shape"):
        pixel_normal(sc.vector([0.0, 0.0, 1.0], unit="m"), "cube")  # type: ignore[arg-type]
