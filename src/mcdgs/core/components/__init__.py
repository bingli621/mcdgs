"""Building blocks of an instrument: one dataclass per McStas component.

Components know nothing about files and import no McStas code. Each provides its McStas
parameters through ``mcstas_parameters`` (its fields plus anything it computes).

A value that follows from a component's own values is a property computed when used. A value that
depends on other components, such as the positions, is left to the instrument.

``COMPONENTS`` maps every component class's name, which is also its McStas component, to the
class. A parameters file names a component's class by it.
"""

from .base import Component
from .choppers import NXdisk_chopper
from .crystal import Single_crystal
from .incoherent import Incoherent
from .masked_source import Polygon_ESS_butterfly
from .misc import Arm, Progress_bar
from .monitors import Monitor_nD
from .optics import Beamstop, Guide_gravity, Pol_bender, Slit
from .phonon import Phonon_simple
from .source import ESS_butterfly

__all__ = [
    "COMPONENTS",
    "Arm",
    "Beamstop",
    "Component",
    "ESS_butterfly",
    "Guide_gravity",
    "Incoherent",
    "Monitor_nD",
    "NXdisk_chopper",
    "Phonon_simple",
    "Pol_bender",
    "Polygon_ESS_butterfly",
    "Progress_bar",
    "Single_crystal",
    "Slit",
]

COMPONENTS: dict[str, type[Component]] = {  # every class exported here but the base
    name: globals()[name] for name in __all__ if name not in ("COMPONENTS", "Component")
}
