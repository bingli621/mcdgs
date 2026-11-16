"""Entities shared by all instruments.

An instrument is a set of components. This package holds no files and no McStas code, and it
imports nothing else from mcdgs. The package layout is described in README.md.

The names below are imported from here, ``from mcdgs.core import Instrument``; the components
from :mod:`.components`.
"""

from .base import INSTRUMENTS, Instrument
from .geometry import place_absolute, rotation
from .params import InstrParam, param
from .utils import bins_from_limits, evaluate, format_value, speed_from_wavelength

__all__ = [
    "INSTRUMENTS",
    "InstrParam",
    "Instrument",
    "bins_from_limits",
    "evaluate",
    "format_value",
    "param",
    "place_absolute",
    "rotation",
    "speed_from_wavelength",
]
