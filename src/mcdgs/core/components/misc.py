"""Components that fit no other module: Progress_bar, Arm."""

from dataclasses import dataclass

from .base import Component


@dataclass(kw_only=True, repr=False)  # repr: Component's
class Progress_bar(Component):
    """McStas ``Progress_bar``, which reports simulation progress.

    Name and instance keywords (``name``, ``at``, ...) are as in :class:`.Component`. The
    McStas parameters left out here (``minutes``, ``flag_save``, ``profile``) keep McStas's
    defaults.

    Parameters
    ----------
    percent:
        Interval between progress reports, in % of the simulation (0 to 100).

    The default is McStas's.

    Examples
    --------
    >>> origin = Progress_bar(name="origin", percent=10)
    >>> origin.mcstas_parameters()
    {'percent': 10}
    """

    percent: float = 10


@dataclass(kw_only=True, repr=False)  # repr: Component's
class Arm(Component):
    """A reference frame, as McStas ``Arm``: it has no McStas parameters.

    An arm whose ``rotated`` uses an instrument parameter acts as a rotation axis.
    Name and instance keywords (``name``, ``at``, ...) are as in :class:`.Component`.

    Examples
    --------
    >>> sample_arm = Arm(name="sample_arm", at=[0, 0, 160], relative="origin")
    >>> sample_arm.mcstas_parameters()
    {}
    """
