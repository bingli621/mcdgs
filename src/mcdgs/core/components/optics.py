"""Neutron optics: guides."""

from dataclasses import dataclass

from .base import Component


@dataclass(kw_only=True, repr=False)  # repr: Component's
class Guide_gravity(Component):
    """A straight or tapered rectangular guide, as McStas ``Guide_gravity``.

    Gravity acts only if ``G`` is set; McStas's default is 0, none.

    Name and instance keywords (``name``, ``at``, ...) are as in :class:`.Component`. The
    McStas parameters left out here (channels, chamfers, sections, rotation, reflectivity file)
    keep McStas's defaults.

    Parameters
    ----------
    w1, h1:
        Width and height at the entrance in m.
    w2, h2:
        Width and height at the exit in m. 0 is the same as at the entrance.
    l:
        Length in m.
    mleft, mright, mtop, mbottom:
        m-value of the left, right, top and bottom mirror. -1 uses ``m``.
    R0:
        Low-angle reflectivity.
    Qc:
        Critical scattering vector in 1/Å.
    alpha:
        Slope of the reflectivity in Å.
    m:
        m-value of all mirrors whose own is -1. 0 is completely absorbing.
    W:
        Width of the supermirror cut-off in 1/Å.
    G:
        Gravitational acceleration in m/s². 0 is no gravity.
    aleft, aright, atop, abottom:
        Slope of the reflectivity (``alpha``) of each mirror. -1 uses ``alpha``.
    wavy:
        Waviness of the mirrors in deg.
    wavy_z, wavy_tb, wavy_lr:
        Waviness in deg along the guide, and across the top/bottom and the left/right mirrors.

    ``w1``, ``h1`` and ``l`` are required, as in McStas; every default is McStas's.

    Examples
    --------
    >>> guide_1 = Guide_gravity(
    ...     w1=0.089955,
    ...     h1=0.029992,
    ...     w2=0.082455,
    ...     h2=0.034937,
    ...     l=0.9994,
    ...     mleft=1.5,
    ...     mright=1.5,
    ...     mtop=4,
    ...     mbottom=4,
    ... )
    >>> guide_1.mcstas_parameters()["m"]
    1.0
    """

    w1: float
    h1: float
    w2: float = 0
    h2: float = 0
    l: float
    mleft: float = -1
    mright: float = -1
    mtop: float = -1
    mbottom: float = -1
    R0: float = 0.995
    Qc: float = 0.0218
    alpha: float = 4.38
    m: float = 1.0
    W: float = 0.003
    G: float = 0
    aleft: float = -1
    aright: float = -1
    atop: float = -1
    abottom: float = -1
    wavy: float = 0
    wavy_z: float = 0
    wavy_tb: float = 0
    wavy_lr: float = 0
