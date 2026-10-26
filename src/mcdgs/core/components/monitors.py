"""Monitors: Monitor_nD."""

from dataclasses import dataclass

from .base import Component


@dataclass(kw_only=True, repr=False)  # repr: Component's
class Monitor_nD(Component):
    """A general monitor of 0D, 1D or 2D signals or event lists, as McStas ``Monitor_nD``.

    What it records, and its shape, are set by ``options`` (see the McStas documentation).
    Name and instance keywords (``name``, ``at``, ...) are as in :class:`.Component`. The
    McStas parameters left out here (user variables ``user0`` to ``user9`` and their names
    ``username0`` to ``username9``) keep McStas's defaults.

    Parameters
    ----------
    xwidth, yheight, zdepth:
        Width, height and thickness of the detector in m. ``xwidth`` is the diameter of a
        disk, sphere, cylinder or banana.
    xmin, xmax, ymin, ymax, zmin, zmax:
        Lower and upper bounds of the opening in m, in place of the width, height and
        thickness.
    bins:
        Number of bins of every variable. 0 uses those in ``options``.
    min, max:
        Range of every variable, in its unit. The defaults use those in ``options``.
    restore_neutron:
        1 leaves the neutron unchanged, as the ``parallel`` option; 0 does not.
    radius:
        Radius of a sphere, cylinder or banana in m. 0 uses ``xwidth``.
    options:
        What to record and the shape, such as ``"lambda limits=[0.5 5] bins=50"``. If ``None``,
        it is unset.
    filename:
        Output file name, in place of the ``file=`` option. If ``None``, McStas names it.
    geometry:
        OFF or PLY file of a detector of any shape. If ``None``, there is none.
    nowritefile:
        1 writes no file; 0 does.
    nexus_bins:
        NeXus output only: -1 stores no bin information, 0 stores it for event lists, 1 for any
        monitor.

    Every default is McStas's; ``None`` is McStas's unset string, ``"NULL"``.

    Examples
    --------
    >>> lambda_monitor = Monitor_nD(
    ...     name="lambda_monitor",
    ...     xwidth=0.1,
    ...     yheight=0.1,
    ...     options="lambda limits=[0.5 5] bins=50",
    ... )
    >>> lambda_monitor.mcstas_parameters()["options"]
    'lambda limits=[0.5 5] bins=50'
    """

    xwidth: float = 0
    yheight: float = 0
    zdepth: float = 0
    xmin: float = 0
    xmax: float = 0
    ymin: float = 0
    ymax: float = 0
    zmin: float = 0
    zmax: float = 0
    bins: int = 0
    min: float = -1e40
    max: float = 1e40
    restore_neutron: int = 0
    radius: float = 0
    options: str | None = None
    filename: str | None = None
    geometry: str | None = None
    nowritefile: int = 0
    nexus_bins: int = 0
