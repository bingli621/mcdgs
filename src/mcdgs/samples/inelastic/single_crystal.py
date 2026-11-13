"""Single crystals that scatter inelastically: lead, from one phonon branch."""

from typing import Any

from mcdgs.core.components import Phonon_simple

__all__ = ["lead"]

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
