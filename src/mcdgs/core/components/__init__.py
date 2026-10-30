"""Building blocks of an instrument: one dataclass per McStas component.

Components know nothing about files and import no McStas code. Each provides its McStas
parameters through ``mcstas_parameters`` (its fields plus anything it computes).

A value that follows from a component's own values is a property computed when used. A value that
depends on other components, such as the positions, is left to the instrument.
"""

from .base import Component
from .crystal import Single_crystal
from .misc import Arm, Progress_bar
from .monitors import Monitor_nD
from .optics import Guide_gravity, Pol_bender
from .source import ESS_butterfly

__all__ = [
    "Arm",
    "Component",
    "ESS_butterfly",
    "Guide_gravity",
    "Monitor_nD",
    "Pol_bender",
    "Progress_bar",
    "Single_crystal",
]
