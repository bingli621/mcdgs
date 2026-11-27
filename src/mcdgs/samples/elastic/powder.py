"""Powders (and other incoherent scatterers) that scatter elastically: Al2O3 and vanadium."""

from typing import Any

from mcdgs.core.components import Incoherent, PowderN

__all__ = ["al2o3", "vanadium"]


def al2o3(**params: Any) -> PowderN:
    """Create an Al2O3 powder (corundum, R-3c) that gives its Debye-Scherrer rings.

    The vanadium cylinder of :func:`vanadium`: 20 mm across and 10 mm high. The reflections
    are the McStas built-in list Al2O3_sapphire.laz, which needs no file. McStas reads the cell
    and the cross-sections from its header: Vc 254.52 Å³, with 6 Al2O3 per cell, and per Al2O3
    a sigma_abs of 0.4625 barns and a sigma_inc of 0.0188 barns. They are left at 0 here, so
    McStas takes them from there: giving them would override the header, and McStas would still
    multiply them by the 6 Al2O3 of the header. The lines are ideal (no width) and each
    neutron scatters once at most.

    No focusing is set: it depends on the instrument. Give the vertical range of the detector
    seen from the sample, e.g. ``al2o3(d_phi=30)`` for a detector ±15 deg about the horizontal
    plane.

    Parameters
    ----------
    **params:
        Further ``PowderN`` parameters, or changes to the defaults (e.g. ``name``, ``at``,
        ``relative``, the shape, ``pack``, ``d_phi``).

    Returns
    -------
    :
        The powder.
    """
    return PowderN(
        **{
            "reflections": "Al2O3_sapphire.laz",
            "radius": 0.01,
            "yheight": 0.01,
            **params,
        }
    )


def vanadium(**params: Any) -> Incoherent:
    """Create a vanadium sample that scatters elastically in every direction.

    A solid cylinder 20 mm across and 10 mm high. The cross-sections are per atom, and the
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
