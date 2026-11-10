"""The chopper cascade of an instrument, without its optics.

:mod:`.builder` builds it as a tof model and gives its disks as scippneutron disk choppers;
:mod:`.acceptance` draws its chopper acceptance diagram, with scippneutron's chopper cascade;
:mod:`.time_distance` draws its time-distance diagram, from the tof model; and :mod:`.reader`
counts the neutrons of a tof run on the bins of a McStas monitor.
"""
