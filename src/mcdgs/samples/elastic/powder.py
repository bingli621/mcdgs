"""Powders (and other incoherent scatterers) that scatter elastically: vanadium."""

from typing import Any

from mcdgs.core.components import Incoherent

__all__ = ["vanadium"]


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
