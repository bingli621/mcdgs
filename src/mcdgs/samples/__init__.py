"""Ready-made samples, for any instrument.

:func:`aluminum` and :func:`la2ni7` return a :class:`.Single_crystal` with its own lattice,
orientation, size and mosaic, and its CIF file (if any) in this folder: they scatter
elastically. :func:`lead` returns a :class:`.Phonon_simple`, which scatters from one phonon
branch. :func:`vanadium` returns an :class:`.Incoherent`, which scatters elastically in every
direction. Placement is up to the caller: ``at``, and ``relative`` to a component of the
instrument, e.g. the goniometer's ``"gonio_phi"``. Keyword arguments set the sample's other
parameters or change the defaults, e.g.
``la2ni7(name="sample", mosaic=60, relative="gonio_phi")``.

``reflections`` names a CIF file relative to this folder (``CIF_DIR`` in
core/components/crystal.py), so a saved parameters file is portable.

To add a sample:

1. Put its CIF file in this folder (e.g. ``NaCl.cif``), or use a McStas built-in reflection list
   (e.g. ``"Al.lau"``), which needs no file.
2. Add a function like :func:`la2ni7`, with ``lattice_params``, ``along_beam``, ``in_plane``,
   ``reflections`` and any other ``Single_crystal`` value it needs (e.g. ``mosaic``), then
   ``**params`` last so a caller can change any of them.
3. Add its name to ``__all__``.
4. Add it to the loop of ``test_the_samples_are_crystals_placed_like_any_component``
   (tests/samples/test_samples.py).
"""

from typing import Any

from mcdgs.core.components import Incoherent, Phonon_simple, Single_crystal

__all__ = ["aluminum", "la2ni7", "lead", "vanadium"]

# a 5 mm cube, single scattering only (much faster)
SMALL_CUBE = {"xwidth": 0.005, "yheight": 0.005, "zdepth": 0.005, "order": 1}


def aluminum(**params: Any) -> Single_crystal:
    """Create an fcc aluminum crystal with cubic axes along the lab axes.

    A 5 mm cube, single scattering only, with a 5' mosaic. Uses the McStas built-in
    reflection list Al.lau.

    Parameters
    ----------
    **params:
        Further ``Single_crystal`` parameters, or changes to the defaults (e.g. ``name``,
        ``at``, ``relative``, ``mosaic``).

    Returns
    -------
    :
        The crystal.
    """
    return Single_crystal(
        **{
            "lattice_params": [4.0495, 4.0495, 4.0495, 90.0, 90.0, 90.0],
            "along_beam": [0, 0, 1],
            "in_plane": [1, 0, 0],
            "reflections": "Al.lau",
            **SMALL_CUBE,
            "mosaic": 5.0,
            **params,
        }
    )


def la2ni7(**params: Any) -> Single_crystal:
    """Create a La2Ni7 crystal (P6_3/mmc) with (001) along the beam.

    The structure is from La2Ni7.cif (COD 1523044). (1.357, 1, 0) lies in the horizontal plane,
    to the left of the beam. A 5 mm cube, single scattering only, with a 30' mosaic.

    Parameters
    ----------
    **params:
        Further ``Single_crystal`` parameters, or changes to the defaults (e.g. ``name``,
        ``at``, ``relative``, ``mosaic``).

    Returns
    -------
    :
        The crystal.
    """
    return Single_crystal(
        **{
            "lattice_params": [5.0556, 5.0556, 24.5980, 90.0, 90.0, 120.0],
            "along_beam": [0, 0, 1],
            "in_plane": [1.357, 1, 0],
            "reflections": "La2Ni7.cif",
            **SMALL_CUBE,
            "mosaic": 30.0,
            **params,
        }
    )


# lead: the velocity of sound along a cubic axis, from the elastic constants and the density
LEAD_DENSITY = 11340.0  # kg/m³
LEAD_C11 = 49.5e9  # Pa, longitudinal along [1 0 0]
LEAD_C44 = 14.9e9  # Pa, transverse along [1 0 0]
MEV_ANGSTROM = 151.93  # m/s, 1 meV Å / hbar: a velocity of sound of 1 meV Å


def lead(**params: Any) -> Phonon_simple:
    """Create an fcc lead crystal that scatters from one phonon branch.

    A cylinder 10 mm across and 10 mm high, at 300 K, with the cubic axes along the lab axes.
    The branch is the longitudinal one along a cubic axis: its velocity of sound is
    ``sqrt(C11 / density)`` = 2089 m/s, or 13.75 meV Å, so it reaches 11.1 meV. For the
    transverse branch, ``sqrt(C44 / density)`` = 1146 m/s, give ``c=7.54``.

    The model has forces between nearest neighbors only, so it follows lead at small ``q`` and
    is too high towards the zone boundary, where lead's longitudinal branch stays below about
    9 meV.

    No focusing window is set: it depends on the instrument. Give the direction of the
    detector from the sample and a window that covers it, e.g.
    ``lead(target_z=1, focus_aw=60, focus_ah=10)`` for one straight ahead.

    Parameters
    ----------
    **params:
        Further ``Phonon_simple`` parameters, or changes to the defaults (e.g. ``name``,
        ``at``, ``relative``, ``T``, ``c``, the focusing window).

    Returns
    -------
    :
        The sample.
    """
    return Phonon_simple(
        **{
            "radius": 0.005,
            "yheight": 0.01,
            "sigma_abs": 0.171,  # barns, at 2200 m/s
            "sigma_inc": 0.0030,  # barns
            "a": 4.9508,  # Å, at room temperature
            "b": 9.405,  # fm
            "M": 207.2,  # amu
            "c": (LEAD_C11 / LEAD_DENSITY) ** 0.5 / MEV_ANGSTROM,  # 13.75 meV Å
            "DW": 1.0,
            "T": 300.0,
            **params,
        }
    )


def vanadium(**params: Any) -> Incoherent:
    """Create a vanadium sample that scatters elastically in every direction.

    A solid cylinder 10 mm across and 10 mm high. The cross-sections are per atom, and the
    volume per atom is half the bcc cell of 3.024 Å. A neutron scatters every 27 mm on average,
    and as many times as it happens to: give ``order=1`` for single scattering only.

    No focusing window is set: it depends on the instrument. Give the direction of the
    detector from the sample and a window that covers it, e.g.
    ``vanadium(target_z=1, focus_aw=60, focus_ah=10)`` for one straight ahead.

    Parameters
    ----------
    **params:
        Further ``Incoherent`` parameters, or changes to the defaults (e.g. ``name``, ``at``,
        ``relative``, the shape, ``pack``, the focusing window).

    Returns
    -------
    :
        The sample.
    """
    return Incoherent(
        **{
            "radius": 0.01,
            "yheight": 0.01,
            "sigma_abs": 5.08,  # barns, at 2200 m/s
            "sigma_inc": 5.08,  # barns
            "Vc": 13.827,  # Å³
            **params,
        }
    )
