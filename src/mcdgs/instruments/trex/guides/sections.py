"""Sections of the T-Rex guide, placed one after the other along the beam axis.

The curved and the straight guides are chains of sections. The sheet's positions are distances
along the beam axis there, so the distance between two entrances is the difference of their
positions. The extraction guides are not a chain: see :mod:`.extraction`.
"""

import math
from dataclasses import dataclass
from itertools import pairwise
from typing import Any

from mcdgs.core.components import Arm, Component, Guide_gravity


@dataclass(frozen=True)
class Section:
    """A guide or a gap of the sheet, with where it lies on the beam axis.

    Parameters
    ----------
    z:
        Start along the beam axis, in m.
    length:
        Length in m.
    component:
        The guide, or the arm at the start of the gap.
    """

    z: float
    length: float
    component: Arm | Guide_gravity

    @property
    def middle(self) -> float:
        """Midpoint along the beam axis, in m."""
        return self.z + self.length / 2


def guide(z: float, l: float, **values: Any) -> Section:
    """Create the section of a guide.

    Parameters
    ----------
    z:
        Start along the beam axis, in m.
    l:
        Length in m.
    values:
        The other keywords of :class:`.Guide_gravity`.

    Returns
    -------
    :
        A new section.
    """
    return Section(z, l, Guide_gravity(l=l, **values))


def gap(z: float, length: float) -> Section:
    """Create the section of a gap.

    Parameters
    ----------
    z:
        Start along the beam axis, in m.
    length:
        Length in m.

    Returns
    -------
    :
        A new section.
    """
    return Section(z, length, Arm())


def arc_angle(distance: float, radius: float) -> float:
    """Compute the turn of the beam axis along an arc.

    Parameters
    ----------
    distance:
        Distance along the arc, in m.
    radius:
        Radius of the arc, in m.

    Returns
    -------
    :
        Rotation about +y, in deg.
    """
    return math.degrees(distance / radius)


def chain(sections: dict[str, Section], radius: float = math.inf) -> dict[str, Component]:
    """Place each section relative to the one before it. The first keeps its own placement.

    Parameters
    ----------
    sections:
        Sections by name, in beam order.
    radius:
        Radius of curvature of the beam axis, in m. The default is a straight axis.

    Returns
    -------
    :
        Components by name, in beam order.
    """
    for (previous_name, previous), (_, section) in pairwise(sections.items()):
        section.component.at = [0, 0, section.z - previous.z]
        section.component.rotated = [0, arc_angle(section.middle - previous.middle, radius), 0]
        section.component.relative = previous_name
    return {name: section.component for name, section in sections.items()}
