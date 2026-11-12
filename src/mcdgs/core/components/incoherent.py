"""Incoherent scatterers: Incoherent, an isotropic elastic or quasielastic sample.

Ready-made incoherent samples are in samples/.
"""

from dataclasses import dataclass

from .base import Component


@dataclass(kw_only=True, repr=False)  # repr: Component's
class Incoherent(Component):
    """An incoherent scatterer such as vanadium, as McStas ``Incoherent``.

    It scatters in every direction with the same probability, whatever the orientation of the
    sample, so it suits a powder as well as a solid piece. The scattering is elastic, with a
    fraction ``f_QE`` of it quasielastic, or it changes the energy by ``Etrans``. The beam is
    attenuated by the absorption, which goes as 1 / velocity, and by the scattering.

    The shape is a box (``xwidth``, ``yheight``, ``zdepth``), a cylinder (``radius``,
    ``yheight``), a sphere (``radius`` alone) or the one of a ``geometry`` file. A
    ``thickness`` makes the box, the cylinder or the sphere hollow.

    The number of atoms per volume is ``pack / Vc``, so the cross-sections and ``Vc`` are for
    the same unit: one atom, or one cell.

    Name and instance keywords (``name``, ``at``, ...) are as in :class:`.Component`.

    Parameters
    ----------
    geometry:
        Name of an OFF or PLY file with the shape. McStas must then be compiled with
        ``-DUSE_OFF``. If ``None``, the shape follows from the sizes.
    radius:
        Outer radius of the cylinder or the sphere, in m.
    xwidth:
        Width of the box, or of the bounding box of ``geometry``, in m.
    yheight:
        Height of the box or the cylinder, or of the bounding box of ``geometry``, in m.
    zdepth:
        Depth of the box, or of the bounding box of ``geometry``, in m.
    thickness:
        Thickness of the wall of a hollow sample, in m. 0 is a solid one.
    target_x, target_y, target_z:
        Position of the target to scatter towards, in m, in the frame of the component.
    target_index:
        Relative index of the component to scatter towards instead: 1 is the next one. McStas
        takes 1 if neither it nor a position is given.
    focus_r:
        Radius of the focusing window, a disk at the target, in m. 0 is every direction.
    focus_xw, focus_yh:
        Width and height of a rectangular focusing window at the target instead, in m.
    focus_aw, focus_ah:
        Width and height of a rectangular focusing window instead, as angles seen from the
        sample, in deg.
    pack:
        Packing factor: the fraction of the volume that is filled, below 1 for a powder.
    p_interact:
        Fraction of the neutrons that is made to scatter; the rest are transmitted. The weights
        make up for it. A value outside (0, 1] leaves it to the cross-sections.
    f_QE:
        Fraction of the scattering that is quasielastic; the rest is elastic.
    gamma:
        Half width at half maximum of the Lorentzian of the quasielastic scattering, in meV.
    Etrans:
        Energy the neutron loses, in meV, the same for every neutron. It is used only with a
        ``deltaE``.
    deltaE:
        Half width of the uniform spread of the energy around ``Etrans``, in meV. 0 is elastic
        and quasielastic scattering.
    sigma_abs:
        Absorption cross-section at 2200 m/s, in barns.
    sigma_inc:
        Incoherent scattering cross-section, in barns.
    Vc:
        Volume that ``sigma_abs`` and ``sigma_inc`` are for, in Å³.
    concentric:
        1 if the sample is hollow and holds other components. It is then given a second time,
        after them. For a box, a cylinder or a sphere only.
    order:
        Maximum scattering order. 0 is unlimited; 1 is single scattering only.

    Notes
    -----
    A size is required: McStas stops without one. Without a focusing window the neutrons are
    scattered in every direction, so few reach a detector: set one that covers it. Every
    default is McStas's, so the cross-sections and ``Vc`` are those of vanadium, per atom.

    Examples
    --------
    >>> sample = Incoherent(name="sample", radius=0.005, yheight=0.01, focus_aw=10, focus_ah=10)
    >>> round(sample.scattering_length, 4)  # m
    0.0272
    """

    geometry: str | None = None
    radius: float = 0
    xwidth: float = 0
    yheight: float = 0
    zdepth: float = 0
    thickness: float = 0
    target_x: float = 0
    target_y: float = 0
    target_z: float = 0
    target_index: int = 0
    focus_r: float = 0
    focus_xw: float = 0
    focus_yh: float = 0
    focus_aw: float = 0
    focus_ah: float = 0
    pack: float = 1
    p_interact: float = 1
    f_QE: float = 0
    gamma: float = 0
    Etrans: float = 0
    deltaE: float = 0
    sigma_abs: float = 5.08
    sigma_inc: float = 5.08
    Vc: float = 13.827
    concentric: int = 0
    order: int = 0

    @property
    def scattering_length(self) -> float:
        """Mean free path between two scatterings, in m: ``Vc / (pack sigma_inc)``."""
        return self.Vc / (self.pack * self.sigma_inc * 100)  # 1 barn / Å³ = 100 / m
