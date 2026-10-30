"""Neutron optics: guides and benders."""

from dataclasses import dataclass, field

from .base import Component


@dataclass(kw_only=True, repr=False)  # repr: Component's
class Guide_gravity(Component):
    """A straight or tapered rectangular guide, as McStas ``Guide_gravity``.

    Name and instance keywords (``name``, ``at``, ...) are as in :class:`.Component`. The
    McStas parameters left out here (channels ``nslit``, ``d``, ``nhslit``; chamfers
    ``chamfers``, ``chamfers_z``, ``chamfers_lr``, ``chamfers_tb``; sections ``nelements``;
    rotation ``nu``, ``phase``; reflectivity file ``reflect``) keep McStas's defaults.

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
        Slope of the reflectivity in Å, of every mirror whose own is -1.
    m:
        m-value of every mirror whose own is -1. 0 is completely absorbing.
    W:
        Width of the supermirror cut-off in 1/Å.
    G:
        Gravitational acceleration in m/s². 0 is no gravity.
    aleft, aright, atop, abottom:
        Slope of the reflectivity of the left, right, top and bottom mirror. -1 uses ``alpha``.
    wavy:
        Waviness of the mirrors in deg.
    wavy_z, wavy_tb, wavy_lr:
        Waviness in deg along the guide, and across the top/bottom and the left/right mirrors.

    ``w1``, ``h1`` and ``l`` are required, as in McStas. Every default is McStas's.

    Examples
    --------
    >>> guide_1 = Guide_gravity(
    ...     name="guide_1",
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


def _default_mirror() -> list[float]:
    """Return McStas's default mirror parameters of ``Pol_bender``: R0, Qc, alpha, m, W."""
    return [0.99, 0.0219, 6.07, 2.0, 0.003]


@dataclass(kw_only=True, repr=False)  # repr: Component's
class Pol_bender(Component):
    """A curved multi-channel guide that can polarise, as McStas ``Pol_bender``.

    Name and instance keywords (``name``, ``at``, ...) are as in :class:`.Component`. The
    McStas parameters left out here (reflectivity files ``rTopUpData`` ... ``rRightDownData``;
    ``debug``) keep McStas's defaults.

    Parameters
    ----------
    xwidth, yheight:
        Width and height at the entrance in m.
    length:
        Length along the centre in m.
    radius:
        Radius of curvature in m. Positive curves left, negative curves right.
    G:
        Gravitational acceleration in m/s². 0 is no gravity.
    nslit:
        Number of slits.
    d:
        Width of the spacers between the slits in m.
    endFlat:
        If positive, the entrance and exit planes are parallel.
    rTopUpPar, rTopDownPar, rBotUpPar, rBotDownPar:
        Reflectivity of the top and bottom mirrors for spin up and down, as
        ``[R0, Qc, alpha, m, W]``.
    rLeftUpPar, rLeftDownPar, rRightUpPar, rRightDownPar:
        Reflectivity of the left and right mirrors for spin up and down, as
        ``[R0, Qc, alpha, m, W]``.
    drawOption:
        Drawing detail: 1 is fine, 2 normal, 3 rough.

    ``xwidth``, ``yheight``, ``length`` and ``radius`` are required, as in McStas. Every default
    is McStas's. Identical up and down values make the bender spin-independent.

    Examples
    --------
    >>> bender = Pol_bender(name="bender", xwidth=0.06, yheight=0.04, length=0.05, radius=7.2)
    >>> bender.mcstas_parameters()["rTopUpPar"]
    [0.99, 0.0219, 6.07, 2.0, 0.003]
    """

    xwidth: float
    yheight: float
    length: float
    radius: float
    G: float = 9.8
    nslit: int = 1
    d: float = 0.0
    endFlat: int = 0
    rTopUpPar: list[float] = field(default_factory=_default_mirror)
    rTopDownPar: list[float] = field(default_factory=_default_mirror)
    rBotUpPar: list[float] = field(default_factory=_default_mirror)
    rBotDownPar: list[float] = field(default_factory=_default_mirror)
    rLeftUpPar: list[float] = field(default_factory=_default_mirror)
    rLeftDownPar: list[float] = field(default_factory=_default_mirror)
    rRightUpPar: list[float] = field(default_factory=_default_mirror)
    rRightDownPar: list[float] = field(default_factory=_default_mirror)
    drawOption: int = 1
