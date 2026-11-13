"""Single crystals that scatter elastically: aluminum and La2Ni7."""

from typing import Any

from mcdgs.core.components import Single_crystal

__all__ = ["aluminum", "la2ni7"]


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
            "xwidth": 0.005,
            "yheight": 0.005,
            "zdepth": 0.005,
            "order": 1,
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
            "reflections": "elastic/La2Ni7.cif",
            "xwidth": 0.005,
            "yheight": 0.005,
            "zdepth": 0.005,
            "order": 1,
            "mosaic": 30.0,
            **params,
        }
    )
