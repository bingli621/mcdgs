"""The source of T-Rex."""

from mcdgs.core.components import (
    ESS_butterfly,
)
from mcdgs.core.components.source import ESS_SOURCE_FREQUENCY, ESS_SOURCE_LENGTH, L_MAX, L_MIN
from mcdgs.core.utils import speed_from_wavelength

SAMPLE_Z = 163.8  # m


def wavelength_band(lambda0: float) -> tuple[float, float]:
    """Return the wavelength band the source simulates around ``lambda0``.

    The band is the wavelengths that reach the sample in one source period plus the pulse
    length, 1.810 Å, centered on ``lambda0`` and clipped to ``L_MIN`` and ``L_MAX``.

    Parameters
    ----------
    lambda0:
        Central wavelength, in Å.

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
    # h / m_n = 3956.034 m·Å/s, Planck's constant over the neutron mass: a neutron of
    # wavelength λ in Å flies at h / (m_n λ) m/s, so the one of 1 Å at 3956.034 m/s
    h_over_m_n = speed_from_wavelength(1.0)  # m·Å/s
    # t = L λ m_n / h, so Δλ = (h / m_n) (1 / f + pulse length) / L. With L = 163.8 m that
    # is 1.810 Å, so 0.905 Å on each side of lambda0
    half_band = 0.5 * h_over_m_n * (1 / ESS_SOURCE_FREQUENCY + ESS_SOURCE_LENGTH) / SAMPLE_Z  # Å

    l_min = max(L_MIN, lambda0 - half_band)
    l_max = min(L_MAX, lambda0 + half_band)
    if l_min >= l_max:
        raise ValueError(f"lambda0 = {lambda0} Å leaves no band between {L_MIN} and {L_MAX} Å")
    return (l_min, l_max)


def source(l_min: float = L_MIN, l_max: float = L_MAX) -> ESS_butterfly:
    """Create the source component for T-Rex.

    Returns
    -------
    :
        A new source.
    """
    return ESS_butterfly(
        sector="W",
        beamline=7,
        yheight=0.03,
        cold_frac=0.3,
        # the focusing rectangle is the entrance of guide 1, in its plane: seen from every point
        # of the moderator it then covers the whole opening, so the flux into guide 1 is right
        dist=1.895927,  # m, cell B13, the entrance of guide 1
        # centered on the axis, it reaches the fixed wall of guide 1 at x = -0.047 m
        # on both sides, so it covers the entrance 2.0 mm off the axis
        focus_xw=0.094,
        focus_yh=0.029992,  # m, cell F13, the height of the entrance of guide 1
        Lmin=l_min,
        Lmax=l_max,
        n_pulses=1,
        acc_power=2,
    )
