"""Phonon samples: Phonon_simple, one acoustic branch of an fcc crystal.

Ready-made phonon samples are in samples/.
"""

from dataclasses import dataclass

from .base import Component


@dataclass(kw_only=True, repr=False)  # repr: Component's
class Phonon_simple(Component):
    """One acoustic phonon branch of an fcc crystal, as McStas ``Phonon_simple``.

    A cylinder of an fcc crystal with one atom per cell and forces between nearest neighbors
    only. Its one branch has the energy ``c / a * sqrt(12 - J(q))``, with ``J(q)`` the sum of
    ``cos(q · d)`` over the twelve nearest neighbors ``d``: ``c |q|`` for a small ``q``, and at
    most ``4 c / a``, at the zone boundary along a cubic axis. It scatters every neutron once,
    into the focusing window, and weights it by the one-phonon cross-section of Squires,
    ch. 3, with the Bose factor of ``T`` for the creation or the annihilation of a phonon.

    The cubic axes of the crystal are the axes of the component: ``rotated`` orients it.

    It leaves out the polarization of the phonons, so every ``q`` scatters whatever its
    direction to the branch; Bragg and incoherent scattering, so there is no elastic line; and
    multiple scattering. ``sigma_abs`` and ``sigma_inc`` only attenuate the beam.

    Name and instance keywords (``name``, ``at``, ...) are as in :class:`.Component`.

    Parameters
    ----------
    radius:
        Radius of the cylinder, in m.
    yheight:
        Height of the cylinder, in m.
    sigma_abs:
        Absorption cross-section per atom at 2200 m/s, in barns.
    sigma_inc:
        Incoherent scattering cross-section per atom, in barns.
    a:
        Lattice constant of the fcc cell, in Å.
    b:
        Coherent scattering length, in fm.
    M:
        Atomic mass, in amu.
    c:
        Velocity of sound, in meV Å: the slope of the branch at small ``q``. 1 meV Å is
        151.9 m/s.
    DW:
        Debye-Waller factor, the same for every ``q``.
    T:
        Temperature, in K.
    target_x, target_y, target_z:
        Position of the target to scatter towards, in m, in the frame of the component.
    target_index:
        Relative index of the component to scatter towards instead: 1 is the next one. McStas
        takes 1 if neither it nor a position is given.
    focus_r:
        Radius of the focusing window, a disk at the target, in m.
    focus_xw, focus_yh:
        Width and height of a rectangular focusing window at the target instead, in m.
    focus_aw, focus_ah:
        Width and height of a rectangular focusing window instead, as angles seen from the
        sample, in deg.
    gap:
        Energy gap of the branch at ``q`` = 0, in meV. It is not physical for a phonon.
    e_steps_low, e_steps_high:
        Number of steps of the search for the final energies on the branch, below and above
        the elastic line.

    Notes
    -----
    The size, the cross-sections, ``a``, ``b``, ``M``, ``c``, ``DW`` and ``T`` are required, as
    in McStas. Without a focusing window the neutrons are scattered in every direction, so few
    reach a detector: set one that covers it. Every default is McStas's.

    Examples
    --------
    >>> sample = Phonon_simple(
    ...     name="sample",
    ...     radius=0.005,
    ...     yheight=0.01,
    ...     sigma_abs=0.171,
    ...     sigma_inc=0.003,
    ...     a=4.9508,
    ...     b=9.405,
    ...     M=207.2,
    ...     c=13.75,
    ...     DW=1,
    ...     T=300,
    ...     focus_aw=10,
    ...     focus_ah=10,
    ... )
    >>> round(sample.branch_maximum, 2)  # meV
    11.11
    """

    radius: float
    yheight: float
    sigma_abs: float
    sigma_inc: float
    a: float
    b: float
    M: float
    c: float
    DW: float
    T: float
    target_x: float = 0
    target_y: float = 0
    target_z: float = 0
    target_index: int = 0
    focus_r: float = 0
    focus_xw: float = 0
    focus_yh: float = 0
    focus_aw: float = 0
    focus_ah: float = 0
    gap: float = 0
    e_steps_low: int = 50
    e_steps_high: int = 50

    @property
    def branch_maximum(self) -> float:
        """Highest energy of the branch, in meV: ``4 c / a``, without a ``gap``."""
        return 4 * self.c / self.a
