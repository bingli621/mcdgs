"""Samples of a tabulated S(q, w): Isotropic_Sqw, of the norm of q only.

Ready-made samples are in samples/.
"""

from dataclasses import dataclass

from .base import Component


@dataclass(kw_only=True, repr=False)  # repr: Component's
class Isotropic_Sqw(Component):
    """A sample of a tabulated S(|q|, w), as McStas ``Isotropic_Sqw``.

    It scatters from the dynamic structure factor of a file, coherent ``Sqw_coh`` and incoherent
    ``Sqw_inc``, elastically and inelastically, with multiple scattering and absorption. It
    knows only the norm of q, not the vector, so it suits a liquid, a glass or a powder: the
    S(|q|, w) of a powder is that of its crystal averaged over the directions of q.

    A file has three blocks of numbers, separated by comment lines that start with ``#``: the
    values of |q|, in 1/Å, then those of w, the energy the neutron loses, in meV, then S(q, w),
    one row per value of |q|. Both axes must be evenly spaced. Comment lines may set
    ``sigma_abs``, ``sigma_coh``, ``sigma_inc``, ``Temperature``, ``V_rho``, ``density``,
    ``weight`` and ``classical``.

    The shape is a box (``xwidth``, ``yheight``, ``zdepth``), a cylinder (``radius``,
    ``yheight``), a sphere (``radius`` alone) or the one of a ``geometry`` file. A
    ``thickness`` makes the box, the cylinder or the sphere hollow.

    Name and instance keywords (``name``, ``at``, ...) are as in :class:`.Component`. The McStas
    parameters left out here, those of a powder file (``powder_format``, ``powder_Dd``,
    ``powder_DW``, ``powder_Vc``, ``powder_barns``) and ``quantum_correction``, keep McStas's
    defaults.

    Parameters
    ----------
    Sqw_coh:
        Name of the file of the coherent S(q, w). If ``None``, there is no coherent scattering.
    Sqw_inc:
        Name of the file of the incoherent S(q, w). If ``None``, the incoherent scattering is
        elastic and the same in every direction, as from vanadium.
    geometry:
        Name of an OFF, PLY or STL file with the shape. If ``None``, the shape follows from the
        sizes.
    radius:
        Outer radius of the cylinder or the sphere, in m.
    thickness:
        Thickness of the wall of a hollow sample, in m. 0 is a solid one.
    xwidth, yheight, zdepth:
        Size of the box, in m; ``yheight`` is also the height of the cylinder.
    threshold:
        Values of S(q, w) below it are left out.
    order:
        Maximum scattering order. 0 is unlimited; 1 is single scattering only.
    T:
        Temperature of the detailed balance, in K. -1 leaves S(q, w) as it is; -2 guesses the
        temperature from the file, which must then hold both signs of w.
    verbose:
        0 silent, 1 normal, 2 verbose, 3 debug.
    d_phi:
        Height of the scattering, as an angle about the horizontal plane, in deg, with single
        scattering (``order=1``) only. 0 is every direction.
    concentric:
        1 if the sample is hollow and holds other components. It is then given a second time,
        after them. For a box, a cylinder or a sphere only.
    rho:
        Number of scattering units per volume, in 1/Å³, that the cross-sections are for.
    sigma_abs:
        Absorption cross-section at 2200 m/s, in barns. -1 leaves out the absorption.
    sigma_coh:
        Coherent scattering cross-section, in barns. -1 leaves out the coherent scattering.
    sigma_inc:
        Incoherent scattering cross-section, in barns. -1 leaves out the incoherent scattering.
    classical:
        1 if S(q, w) is classical, symmetric in w, so that McStas multiplies it by the detailed
        balance; 0 if it holds the Bose factor already, as measured; -1 guesses from the file.
    density:
        Density of the material, in g/cm³: with ``weight``, instead of ``rho``.
    weight:
        Weight of a scattering unit, in g/mol.
    p_interact:
        Fraction of the neutrons that is made to scatter; the weights make up for it. -1 leaves
        it to the cross-sections.
    norm:
        -1 scales S(q, w) so that S(q) tends to 1 at large q; 0 or 1 takes it as it is; above 0,
        multiplies it.

    Notes
    -----
    Every default is McStas's: give a size, a file or the cross-sections of an incoherent
    scatterer, and ``rho`` or ``density`` and ``weight``. To take S(q, w) as it is, give
    ``T=-1`` and ``norm=1``.

    Examples
    --------
    >>> sample = Isotropic_Sqw(
    ...     name="sample", Sqw_coh="magnon.sqw", radius=0.01, yheight=0.01, T=-1, norm=1
    ... )
    >>> sample.mcstas_parameters()["Sqw_inc"] is None  # NULL in McStas
    True
    """

    Sqw_coh: str | None = None
    Sqw_inc: str | None = None
    geometry: str | None = None
    radius: float = 0
    thickness: float = 0
    xwidth: float = 0
    yheight: float = 0
    zdepth: float = 0
    threshold: float = 1e-20
    order: int = 0
    T: float = 0
    verbose: int = 1
    d_phi: float = 0
    concentric: int = 0
    rho: float = 0
    sigma_abs: float = 0
    sigma_coh: float = 0
    sigma_inc: float = 0
    classical: float = -1
    density: float = 0
    weight: float = 0
    p_interact: float = -1
    norm: float = -1
