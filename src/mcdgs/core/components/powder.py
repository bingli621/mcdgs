"""Powders: PowderN, the Debye-Scherrer rings of a list of reflections.

Ready-made powders are in samples/.
"""

from dataclasses import dataclass
from typing import Any

from .base import Component
from .crystal import CIF_DIR


@dataclass(kw_only=True, repr=False)  # repr: Component's
class PowderN(Component):
    """A powder, as McStas ``PowderN``.

    Each neutron is transmitted, scattered incoherently or Bragg scattered on the cone of one
    reflection of ``reflections``, once at most: there is no multiple scattering. The beam is
    attenuated by the absorption, which goes as 1 / velocity, and by the scattering.

    The shape is a box (``xwidth``, ``yheight``, ``zdepth``), a cylinder (``radius``,
    ``yheight``), a sphere (``radius`` alone) or the one of a ``geometry`` file. A
    ``thickness`` makes the box, the cylinder or the sphere hollow.

    ``Vc``, ``sigma_abs``, ``sigma_inc``, ``DW`` and ``delta_d_d`` left at 0 are read from the
    header of ``reflections`` if it has them, as a LAZ file of McStas does (e.g.
    ``"Al2O3_sapphire.laz"``): a value given here overrides the header. The cross-sections and
    ``Vc`` are per unit cell.

    Name and instance keywords (``name``, ``at``, ...) are as in :class:`.Component`. The McStas
    parameter left out here, ``format``, keeps McStas's default: the columns of ``reflections``
    are those its header names, or those of a LAZ file.

    Parameters
    ----------
    reflections:
        A reflection list (LAZ, LAU, FullProf, ShelX), an NCrystal cfg-string (with ``.ncmat``
        or ``::``), or a CIF file relative to samples/. A name alone is looked for in the data
        folder of McStas. If ``None``, only the incoherent scattering is left.
    geometry:
        Name of an OFF, PLY or STL file with the shape. If ``None``, the shape follows from the
        sizes.
    radius:
        Outer radius of the cylinder or the sphere, in m.
    yheight:
        Height of the box or the cylinder, or of the bounding box of ``geometry``, in m.
    xwidth:
        Width of the box, or of the bounding box of ``geometry``, in m.
    zdepth:
        Depth of the box, or of the bounding box of ``geometry``, in m.
    thickness:
        Thickness of the wall of a hollow sample, in m. 0 is a solid one; a negative value puts
        the hollow outside the box or the cylinder.
    pack:
        Packing factor: the fraction of the volume that is filled.
    Vc:
        Volume of the unit cell, in Å³. 0 reads it from ``reflections``, or computes it from
        ``density``, ``weight`` and ``nb_atoms``.
    sigma_abs:
        Absorption cross-section per unit cell at 2200 m/s, in barns. A negative value leaves
        out the absorption.
    sigma_inc:
        Incoherent scattering cross-section per unit cell, in barns. A negative value leaves
        out the incoherent scattering.
    delta_d_d:
        Relative width Δd/d of every line, unless ``reflections`` has a column of it. 0 is an
        ideal powder.
    p_inc:
        Fraction of the neutrons that is scattered incoherently. The weights make up for it.
    p_transmit:
        Fraction of the neutrons that is only attenuated, not scattered.
    DW:
        Debye-Waller factor of every line, unless ``reflections`` has a column of it. 0 reads it
        from ``reflections``; give 1 if it is in the structure factors already.
    nb_atoms:
        Number of formula units per unit cell: the ratio of the cross-sections per cell to
        those of the chemical formula.
    d_omega:
        Horizontal range to focus the incoherent scattering on, in deg. 0 is every direction.
    d_phi:
        Vertical range of the Debye-Scherrer cones to focus on, in deg, such as the height of
        the detector seen from the sample. 0 is the whole cone.
    tth_sign:
        Side of the Bragg scattering, the sign of the scattering angle. 0 is both sides at
        random. Used only with a ``d_phi``.
    p_interact:
        Fraction of the neutrons that is Bragg scattered. The weights make up for it.
    concentric:
        1 if the sample is hollow and holds other components. It is then given a second time,
        after them. For a box, a cylinder or a sphere only.
    density:
        Density, in g/cm³, to compute ``Vc`` from.
    weight:
        Molar mass of the chemical formula, in g/mol, to compute ``Vc`` from.
    barns:
        1 if the structure factors of ``reflections`` are in barns (LAZ), 0 if in fm² (LAU).
    Strain:
        Relative shift Δd/d of every line, in ppm, unless ``reflections`` has a column of it.
    focus_flip:
        0 if ``d_phi`` is measured from the horizontal x-z plane, 1 if from the vertical y-z
        plane.
    target_index:
        Relative index of the component to focus the incoherent scattering on: 1 is the next
        one. 0 is none.

    Notes
    -----
    A size is required: McStas stops without one. A ``d_phi`` that covers the detector spends
    the neutrons where it is: without it, most go up or down along the cones. Every default is
    McStas's.

    Examples
    --------
    >>> sample = PowderN(
    ...     name="sample", reflections="Al2O3_sapphire.laz", radius=0.005, yheight=0.01, d_phi=30
    ... )
    >>> sample.mcstas_parameters()["reflections"]  # found by McStas in its data folder
    'Al2O3_sapphire.laz'
    """

    reflections: str | None = None
    geometry: str | None = None
    radius: float = 0
    yheight: float = 0
    xwidth: float = 0
    zdepth: float = 0
    thickness: float = 0
    pack: float = 1
    Vc: float = 0
    sigma_abs: float = 0
    sigma_inc: float = 0
    delta_d_d: float = 0
    p_inc: float = 0.1
    p_transmit: float = 0.1
    DW: float = 0
    nb_atoms: float = 1
    d_omega: float = 0
    d_phi: float = 0
    tth_sign: int = 0
    p_interact: float = 0.8
    concentric: int = 0
    density: float = 0
    weight: float = 0
    barns: int = 1
    Strain: float = 0
    focus_flip: int = 0
    target_index: int = 0

    def mcstas_parameters(self) -> dict[str, Any]:
        """Return the McStas parameters.

        Returns
        -------
        :
            The fields, with ``reflections`` given as a full path if it names a CIF file in
            samples/, else unchanged.
        """
        if self.reflections is None or not (cif := CIF_DIR / self.reflections).is_file():
            return super().mcstas_parameters()
        return super().mcstas_parameters() | {"reflections": str(cif)}
