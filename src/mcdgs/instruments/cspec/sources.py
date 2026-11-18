"""The source of CSPEC: :func:`source`, or :func:`masked_source`, masked by the choppers.

Both have the same values, from :func:`_source`: they focus on the entrance of the guide system,
the entrance of W03-01-011, with the same band.
"""

from typing import Any

from mcdgs.core import speed_from_wavelength
from mcdgs.core.components import ESS_butterfly, Polygon_ESS_butterfly
from mcdgs.core.components.masked_source import ChopperTrainParameters
from mcdgs.core.components.source import ESS_SOURCE_FREQUENCY, ESS_SOURCE_LENGTH, L_MAX, L_MIN

from .guides.bender import (
    GUIDE_01_011_IN_HEIGHT,
    GUIDE_01_011_IN_WIDTH,
    GUIDE_START_X,
    GUIDE_START_Z,
)
from .guides.frame import MM_PER_M
from .guides.straight import SAMPLE_PATH


def wavelength_band(lambda0: float | None) -> tuple[float, float]:
    """Return the wavelength band the source simulates around ``lambda0``.

    The band is the wavelengths that reach the sample in one source period plus the pulse
    length, 1.853 Å, centered on ``lambda0`` and clipped to ``L_MIN`` and ``L_MAX``. Without a
    ``lambda0`` it is not limited: the whole of ``L_MIN`` to ``L_MAX``.

    Parameters
    ----------
    lambda0:
        Central wavelength, in Å, or ``None`` for the white beam.

    Returns
    -------
    :
        ``Lmin`` and ``Lmax`` of the source, in Å.

    Raises
    ------
    ValueError
        If the clipped band is empty: ``lambda0`` is half the band or more outside ``L_MIN``
        to ``L_MAX``.
    """
    if lambda0 is None:
        return (L_MIN, L_MAX)
    # h / m_n = 3956.034 m·Å/s, Planck's constant over the neutron mass: a neutron of
    # wavelength λ in Å flies at h / (m_n λ) m/s, so the one of 1 Å at 3956.034 m/s
    h_over_m_n = speed_from_wavelength(1.0)  # m·Å/s
    # t = L λ m_n / h, so Δλ = (h / m_n) (1 / f + pulse length) / L. With L = 160.004 m that
    # is 1.853 Å, so 0.927 Å on each side of lambda0
    half_band = 0.5 * h_over_m_n * (1 / ESS_SOURCE_FREQUENCY + ESS_SOURCE_LENGTH) / SAMPLE_PATH

    l_min = max(L_MIN, lambda0 - half_band)
    l_max = min(L_MAX, lambda0 + half_band)
    if l_min >= l_max:
        raise ValueError(f"lambda0 = {lambda0} Å leaves no band between {L_MIN} and {L_MAX} Å")
    return (l_min, l_max)


# extra flight path a neutron may have, as a fraction of each disk's path: the guides are not a
# straight line. FIXME, not in CSPEC_nl-specs_1.E.xlsx: Polygon_ESS_butterfly's "realistic figure"
PATH_SPREAD_FRACTION = 1e-4


def _source(l_min: float, l_max: float) -> dict[str, Any]:
    """Return the values both sources share.

    Parameters
    ----------
    l_min, l_max:
        Wavelength band, in Å.

    Returns
    -------
    :
        The McStas parameters of ``ESS_butterfly``, by name.
    """
    return {
        "sector": "W",
        "beamline": 3,
        "yheight": 0.03,
        # FIXME, not in CSPEC_nl-specs_1.E.xlsx
        "cold_frac": 0.9,
        # the focusing rectangle is in the plane of the entrance of the guide system: seen from
        # every point of the moderator it then covers the whole opening, so the flux into the
        # guide is right
        "dist": GUIDE_START_Z,
        # centered on the axis, it reaches the far wall of the entrance, 16.3 mm off the axis,
        # on both sides
        "focus_xw": 2 * (GUIDE_START_X + GUIDE_01_011_IN_WIDTH / MM_PER_M / 2),
        "focus_yh": GUIDE_01_011_IN_HEIGHT / MM_PER_M,
        "Lmin": l_min,
        "Lmax": l_max,
        "n_pulses": 1,
        "acc_power": 2,
    }


def source(l_min: float = L_MIN, l_max: float = L_MAX) -> ESS_butterfly:
    """Create the source component for CSPEC.

    Parameters
    ----------
    l_min, l_max:
        Wavelength band, in Å.

    Returns
    -------
    :
        A new source.
    """
    return ESS_butterfly(**_source(l_min, l_max))


def masked_source(
    choppers: list[ChopperTrainParameters], l_min: float = L_MIN, l_max: float = L_MAX
) -> Polygon_ESS_butterfly:
    """Create the source of CSPEC masked by its choppers.

    It is :func:`source` emitting only where the disks of ``choppers`` let neutrons through, with
    :data:`PATH_SPREAD_FRACTION`. Every ray then reaches the choppers in a window they are open
    in, with its weight multiplied by the acceptance.

    Parameters
    ----------
    choppers:
        The chopper train, as :func:`.chopper_train` describes the choppers of the instrument
        with their flight paths.
    l_min, l_max:
        Wavelength band, in Å.

    Returns
    -------
    :
        A new source.
    """
    return Polygon_ESS_butterfly(
        **_source(l_min, l_max),
        choppers=choppers,
        path_spread_fraction=PATH_SPREAD_FRACTION,
    )
