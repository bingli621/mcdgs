"""The monitors of T-Rex.

Each monitor leaves the neutrons unchanged. The source monitor is the entrance of guide 1;
the beam monitor in the middle of the bunker gap of row 29, after guide 4, on the arm of guide 4;
the pre-sample monitor between ``slit_mc`` and the sample; the beamstop monitor on the face of
the beamstop, but before the sample in beam order; the sample monitor at ``sample_pos``.

The beamstop is here too, behind its monitor: an absorbing plate on the beam axis just in front
of the banana detector of :mod:`.detectors`, which stops the direct beam there.
"""

from mcdgs.core.components import Beamstop, Monitor_nD
from mcdgs.core.components.monitors import monitor_time_limits
from mcdgs.core.components.source import ESS_SOURCE_LENGTH, L_MAX, L_MIN
from mcdgs.core.utils import bins_from_limits

from .detectors import BANANA_RADIUS
from .guides.bispectral import GUIDE_4_TILT, GUIDE_4A_IN_Z
from .guides.extraction import guides_extraction, z_along_arm
from .guides.straight import ARM_MC_Z, SAMPLE_Z

SOURCE_MONITOR_Z = 1.895927  # m, cell B13, the entrance of guide 1
SOURCE_MONITOR_LAMBDA_RESOLUTION = 0.05  # Å, the width of a wavelength bin
SOURCE_MONITOR_TIME_RESOLUTION = 5e-5  # s, the width of a time bin, 0.05 ms

BEAM_MONITOR_GAP_Z = 5.892  # m, cell B29
BEAM_MONITOR_GAP_LENGTH = 0.02  # m, cell D29
BEAM_MONITOR_Z = BEAM_MONITOR_GAP_Z + BEAM_MONITOR_GAP_LENGTH / 2  # m, mid-gap
BEAM_MONITOR_LAMBDA_RESOLUTION = 0.02  # Å, the width of a wavelength bin
BEAM_MONITOR_TIME_RESOLUTION = 5e-5  # s, the width of a time bin, 0.05 ms

PRE_SAMPLE_MONITOR_Z = 163.5  # m, FIXME, not in Table of Optics v3.2: between slit and sample
PRE_SAMPLE_MONITOR_WIDTH = 0.1  # m, FIXME, not in Table of Optics v3.2: wider than the beam
PRE_SAMPLE_MONITOR_HEIGHT = 0.1  # m, FIXME, not in Table of Optics v3.2: wider than the beam
PRE_SAMPLE_MONITOR_LAMBDA_RESOLUTION = 0.02  # Å, FIXME, as the beam monitor
PRE_SAMPLE_MONITOR_TIME_RESOLUTION = 5e-5  # s, FIXME, as the beam monitor

BEAMSTOP_GAP = 0.001  # m, FIXME, not in Table of Optics v3.2: in front of the banana detector
BEAMSTOP_DISTANCE = BANANA_RADIUS - BEAMSTOP_GAP  # m, from the sample
BEAMSTOP_WIDTH = PRE_SAMPLE_MONITOR_WIDTH  # m, FIXME, as the pre-sample monitor
BEAMSTOP_HEIGHT = PRE_SAMPLE_MONITOR_HEIGHT  # m, FIXME, as the pre-sample monitor
BEAMSTOP_MONITOR_LAMBDA_RESOLUTION = 0.02  # Å, FIXME, as the beam monitor
BEAMSTOP_MONITOR_TIME_RESOLUTION = 5e-5  # s, FIXME, as the beam monitor


def time_lambda_options(
    distance: float,
    time_limits: tuple[float, float] | None,
    lambda_limits: tuple[float, float] | None,
    lambda_resolution: float,
    time_resolution: float,
) -> str:
    """Return the ``options`` of a monitor that records one 2D map of time against wavelength.

    Parameters
    ----------
    distance:
        Flight path from the moderator to the monitor, in m.
    time_limits:
        Start and stop of the time axis, in s. ``None`` takes them from
        :func:`.monitor_time_limits` at ``distance``, for ``lambda_limits`` and a pulse of
        ``ESS_SOURCE_LENGTH``.
    lambda_limits:
        Start and stop of the wavelength axis, in Å. ``None`` takes the widest band the source
        simulates, ``L_MIN`` to ``L_MAX``.
    lambda_resolution:
        Width of a wavelength bin, in Å. The axis has the fewest bins that cover
        ``lambda_limits``, see :func:`.bins_from_limits`.
    time_resolution:
        Width of a time bin, in s, likewise.

    Returns
    -------
    :
        Such as ``"time limits=[0.0016 0.0078] bins=124 lambda limits=[1.095 2.915] bins=91"``.
    """
    if lambda_limits is None:
        lambda_limits = (L_MIN, L_MAX)
    if time_limits is None:
        time_limits = monitor_time_limits(distance, lambda_limits, ESS_SOURCE_LENGTH)
    t_start, t_stop, t_bins = bins_from_limits(time_limits, time_resolution)
    l_start, l_stop, l_bins = bins_from_limits(lambda_limits, lambda_resolution)
    return (  # time in s, lambda in Å: one 2D map
        f"time limits=[{t_start} {t_stop}] bins={t_bins} "
        f"lambda limits=[{l_start} {l_stop}] bins={l_bins}"
    )


def source_monitor(
    lambda_resolution: float = SOURCE_MONITOR_LAMBDA_RESOLUTION,
    time_resolution: float = SOURCE_MONITOR_TIME_RESOLUTION,
) -> Monitor_nD:
    """Create the time and wavelength monitor at the entrance of guide 1.

    It is the entrance of guide 1: its opening, at the place and with the rotation of
    ``arm_nboa``, 2.0 mm off the axis of the source and turned by the taper of the extraction
    guides, so it counts what enters guide 1. It records one 2D map of time against wavelength, and leaves
    the neutrons unchanged. The wavelength axis is the widest band the source simulates,
    ``L_MIN`` to ``L_MAX``; the time axis is :func:`.monitor_time_limits` at
    ``SOURCE_MONITOR_Z``, for that band and a pulse of ``ESS_SOURCE_LENGTH``.

    Parameters
    ----------
    lambda_resolution:
        Width of a wavelength bin, in Å. The axis has the fewest bins that cover the band, see
        :func:`.bins_from_limits`.
    time_resolution:
        Width of a time bin, in s, likewise.

    Returns
    -------
    :
        A new monitor.
    """
    lambda_limits = (L_MIN, L_MAX)
    time_limits = monitor_time_limits(SOURCE_MONITOR_Z, lambda_limits, ESS_SOURCE_LENGTH)
    t_start, t_stop, t_bins = bins_from_limits(time_limits, time_resolution)
    l_start, l_stop, l_bins = bins_from_limits(lambda_limits, lambda_resolution)
    entrance = guides_extraction()["arm_nboa"]  # the entrance of guide 1
    return Monitor_nD(
        xwidth=0.089955,  # m, cell E13, the entrance of guide 1
        yheight=0.029992,  # m, cell F13, the entrance of guide 1
        restore_neutron=1,
        options=(  # time in s, lambda in Å: one 2D map
            f"time limits=[{t_start} {t_stop}] bins={t_bins} "
            f"lambda limits=[{l_start} {l_stop}] bins={l_bins}"
        ),
        filename="source_monitor",  # the name of the monitor, in both models
        at=list(entrance.at),  # the arm is placed absolute, as the source: the same frame
        rotated=list(entrance.rotated),
        relative="source",
    )


def beam_monitor(
    time_limits: tuple[float, float] | None = None,
    lambda_limits: tuple[float, float] | None = None,
    lambda_resolution: float = BEAM_MONITOR_LAMBDA_RESOLUTION,
    time_resolution: float = BEAM_MONITOR_TIME_RESOLUTION,
) -> Monitor_nD:
    """Create the time and wavelength monitor in the middle of the bunker gap of row 29.

    It records the intensity leaving guide 4 as one 2D map of time against wavelength, and
    leaves the neutrons unchanged.

    Parameters
    ----------
    time_limits:
        Start and stop of the time axis, in s. ``None`` takes them from
        :func:`.monitor_time_limits` at ``BEAM_MONITOR_Z``, for ``lambda_limits`` and a pulse of
        ``ESS_SOURCE_LENGTH``.
    lambda_limits:
        Start and stop of the wavelength axis, in Å. ``None`` takes the widest band the source
        simulates, ``L_MIN`` to ``L_MAX``.
    lambda_resolution:
        Width of a wavelength bin, in Å. The axis has the fewest bins that cover
        ``lambda_limits``, see :func:`.bins_from_limits`.
    time_resolution:
        Width of a time bin, in s, likewise.

    Returns
    -------
    :
        A new monitor.
    """
    return Monitor_nD(
        xwidth=0.06,  # FIXME cell E32, the entrance of guide 5
        yheight=0.046260695185495173,  # FIXME cell F32, the entrance of guide 5
        restore_neutron=1,
        options=time_lambda_options(
            BEAM_MONITOR_Z, time_limits, lambda_limits, lambda_resolution, time_resolution
        ),
        filename="beam_monitor",
        at=[0, 0, z_along_arm(BEAM_MONITOR_Z, GUIDE_4A_IN_Z, GUIDE_4_TILT)],
        relative="arm_guide_4",
    )


def pre_sample_monitor(
    time_limits: tuple[float, float] | None = None,
    lambda_limits: tuple[float, float] | None = None,
    lambda_resolution: float = PRE_SAMPLE_MONITOR_LAMBDA_RESOLUTION,
    time_resolution: float = PRE_SAMPLE_MONITOR_TIME_RESOLUTION,
) -> Monitor_nD:
    """Create the time and wavelength monitor between ``slit_mc`` and the sample.

    It records the beam leaving the slit as one 2D map of time against wavelength, and leaves
    the neutrons unchanged. It is ``PRE_SAMPLE_MONITOR_Z`` along the beam, on the arm of the
    monochromatic choppers.

    Parameters
    ----------
    time_limits, lambda_limits, lambda_resolution, time_resolution:
        As in :func:`time_lambda_options`, at ``PRE_SAMPLE_MONITOR_Z``.

    Returns
    -------
    :
        A new monitor.
    """
    return Monitor_nD(
        xwidth=PRE_SAMPLE_MONITOR_WIDTH,
        yheight=PRE_SAMPLE_MONITOR_HEIGHT,
        restore_neutron=1,
        options=time_lambda_options(
            PRE_SAMPLE_MONITOR_Z, time_limits, lambda_limits, lambda_resolution, time_resolution
        ),
        filename="pre_sample_monitor",
        at=[0, 0, PRE_SAMPLE_MONITOR_Z - ARM_MC_Z],
        relative="arm_mc",
    )


def sample_monitor() -> Monitor_nD:
    """Create the wavelength monitor at the sample position.

    It records the beam reaching the sample and leaves the neutrons unchanged. It is at
    ``sample_pos``, before the sample.

    Returns
    -------
    :
        A new monitor.
    """
    return Monitor_nD(
        xwidth=0.1,  # FIXME, not in Table of Optics v3.2: wider than the beam
        yheight=0.1,  # FIXME, not in Table of Optics v3.2: wider than the beam
        restore_neutron=1,
        options="lambda limits=[0.5 5.3] bins=12",  # FIXME, not in Table of Optics v3.2
        filename="sample_monitor",
        relative="sample_pos",
    )


def beamstop_monitor(
    time_limits: tuple[float, float] | None = None,
    lambda_limits: tuple[float, float] | None = None,
    lambda_resolution: float = BEAMSTOP_MONITOR_LAMBDA_RESOLUTION,
    time_resolution: float = BEAMSTOP_MONITOR_TIME_RESOLUTION,
) -> Monitor_nD:
    """Create the time and wavelength monitor on the face of the beamstop.

    It has the size and the place of the beamstop, ``BEAMSTOP_DISTANCE`` after ``sample_pos`` on
    the beam axis, but comes before the sample in beam order, so it records the direct beam as if
    there were no sample: each neutron is taken to its plane, counted, and given back unchanged
    to the sample. It records one 2D map of time against wavelength. ``sample_pos`` comes after
    it, so it hangs off ``arm_mc``, whose axis runs straight to the sample.

    Parameters
    ----------
    time_limits, lambda_limits, lambda_resolution, time_resolution:
        As in :func:`time_lambda_options`, at ``SAMPLE_Z + BEAMSTOP_DISTANCE`` from the
        moderator.

    Returns
    -------
    :
        A new monitor.
    """
    return Monitor_nD(
        xwidth=BEAMSTOP_WIDTH,
        yheight=BEAMSTOP_HEIGHT,
        restore_neutron=1,
        options=time_lambda_options(
            SAMPLE_Z + BEAMSTOP_DISTANCE,
            time_limits,
            lambda_limits,
            lambda_resolution,
            time_resolution,
        ),
        filename="beamstop_monitor",
        at=[0, 0, SAMPLE_Z + BEAMSTOP_DISTANCE - ARM_MC_Z],
        relative="arm_mc",
    )


def beamstop() -> Beamstop:
    """Create the beamstop, on the beam axis ``BEAMSTOP_DISTANCE`` after ``sample_pos``.

    Returns
    -------
    :
        A new beamstop.
    """
    return Beamstop(
        xwidth=BEAMSTOP_WIDTH,
        yheight=BEAMSTOP_HEIGHT,
        at=[0, 0, BEAMSTOP_DISTANCE],
        relative="sample_pos",
    )
