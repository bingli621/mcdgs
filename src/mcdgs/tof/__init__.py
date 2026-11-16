"""The chopper cascade of an instrument, without its optics.

:mod:`.builder` builds it as a tof model and gives its disks as scippneutron disk choppers;
:mod:`.acceptance` draws its chopper acceptance diagram, with scippneutron's chopper cascade;
:mod:`.time_distance` draws its time-distance diagram, from the tof model; and :mod:`.reader`
counts the neutrons of a tof run on the bins of a McStas monitor.

The names below are imported from here, ``from mcdgs.tof import build_tof``.
"""

from .acceptance import acceptance_diagram, chopper_cascade
from .builder import beamline, build_tof, disk_chopper, get_choppers
from .reader import monitor_bins, read_tof_monitor
from .time_distance import time_distance_diagram

__all__ = [
    "acceptance_diagram",
    "beamline",
    "build_tof",
    "chopper_cascade",
    "disk_chopper",
    "get_choppers",
    "monitor_bins",
    "read_tof_monitor",
    "time_distance_diagram",
]
