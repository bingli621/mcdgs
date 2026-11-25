"""The monitors of CSPEC.

Each monitor records one 2D map of time against wavelength and leaves the neutrons unchanged.
In the order of the beam, at the places of the monitors of drcspec:

- ``source_monitor``: the entrance of W03-01-011, where the guide system begins;
- ``beam_monitor``: the middle of the gap after the BBG, row 22, 5.9 m from the moderator;
- ``bwfa_monitor``: the middle of the gap after the BWFA, row 50, which the sheet keeps "for
  monitor etc.", 27.9 m;
- ``psc_monitor``: in the gap of the pulse-shaping choppers, halfway from ``ps2`` to its end,
  105.7 m;
- ``mc_monitor``: in the gap of the monochromatic choppers, halfway from ``mc2`` to its end,
  158.5 m;
- ``sample_monitor``: at ``sample_pos``, before the sample, 160.0 m.

CSPEC_nl-specs_1.E.xlsx has no monitors. Each monitor in a gap is as large as the entrance of
the guide after it, and the bins are those of the monitors of T-Rex.
"""

from mcdgs.core.components import Monitor_nD
from mcdgs.core.components.monitors import time_lambda_options

from .choppers import MC_DISK_SPACING, MC_GAP_BEFORE, PSC_DISK_SPACING, PSC_GAP_BEFORE
from .guides.bender import (
    ARM_GAP_BBG_PATH,
    ARM_GAP_BWFA_PATH,
    GAP_BBG_LENGTH,
    GAP_BWFA_LENGTH,
    GUIDE_01_011_IN_HEIGHT,
    GUIDE_01_011_IN_WIDTH,
    GUIDE_START_Z,
    guides_bender,
)
from .guides.frame import MM_PER_M
from .guides.straight import ARM_MC_PATH, ARM_PSC_PATH, SAMPLE_PATH

SOURCE_MONITOR_LAMBDA_RESOLUTION = 0.05  # Å, the width of a wavelength bin, as for T-Rex
SOURCE_MONITOR_TIME_RESOLUTION = 5e-5  # s, the width of a time bin, 0.05 ms, as for T-Rex

# the monitors in the gaps: FIXME, not in CSPEC_nl-specs_1.E.xlsx
GAP_MONITOR_LAMBDA_RESOLUTION = 0.02  # Å, FIXME, as the beam monitor of T-Rex
GAP_MONITOR_TIME_RESOLUTION = 5e-5  # s, FIXME, as the beam monitor of T-Rex

BEAM_MONITOR_Z = GAP_BBG_LENGTH / 2  # mm, from the start of the gap, its middle
BEAM_MONITOR_PATH = ARM_GAP_BBG_PATH + BEAM_MONITOR_Z / MM_PER_M  # m
BEAM_MONITOR_WIDTH = 70  # mm, FIXME cell X23, the entrance of W03-03-01
BEAM_MONITOR_HEIGHT = 55  # mm, FIXME cell Z23, the entrance of W03-03-01

BWFA_MONITOR_Z = GAP_BWFA_LENGTH / 2  # mm, from the start of the gap, its middle
BWFA_MONITOR_PATH = ARM_GAP_BWFA_PATH + BWFA_MONITOR_Z / MM_PER_M  # m
BWFA_MONITOR_WIDTH = 70  # mm, FIXME cell X51, the entrance of W03-06-011
BWFA_MONITOR_HEIGHT = 55  # mm, FIXME cell Z51, the entrance of W03-06-011

PSC_GAP_LENGTH = 107.31  # mm, cell T134
# halfway from ps2 to the end of the gap. FIXME, not in CSPEC_nl-specs_1.E.xlsx
PSC_MONITOR_Z = (PSC_GAP_BEFORE + PSC_DISK_SPACING + PSC_GAP_LENGTH) / 2  # mm
PSC_MONITOR_PATH = ARM_PSC_PATH + PSC_MONITOR_Z / MM_PER_M  # m
PSC_MONITOR_WIDTH = 114.66  # mm, FIXME cell X135, the entrance of W03-08-01
PSC_MONITOR_HEIGHT = 73.89  # mm, FIXME cell Z135, the entrance of W03-08-01

MC_GAP_LENGTH = 72.97  # mm, cell T188
# halfway from mc2 to the end of the gap. FIXME, not in CSPEC_nl-specs_1.E.xlsx
MC_MONITOR_Z = (MC_GAP_BEFORE + MC_DISK_SPACING + MC_GAP_LENGTH) / 2  # mm
MC_MONITOR_PATH = ARM_MC_PATH + MC_MONITOR_Z / MM_PER_M  # m
MC_MONITOR_WIDTH = 22.58  # mm, FIXME cell X189, the entrance of W03-09-02
MC_MONITOR_HEIGHT = 49  # mm, FIXME cell Z189, the entrance of W03-09-02

SAMPLE_MONITOR_WIDTH = 0.1  # m, FIXME, not in CSPEC_nl-specs_1.E.xlsx: wider than the beam
SAMPLE_MONITOR_HEIGHT = 0.1  # m, FIXME, not in CSPEC_nl-specs_1.E.xlsx: wider than the beam
SAMPLE_MONITOR_LAMBDA_RESOLUTION = 0.02  # Å, FIXME, as the sample monitor of T-Rex
SAMPLE_MONITOR_TIME_RESOLUTION = 5e-5  # s, FIXME, as the sample monitor of T-Rex


def source_monitor(
    lambda_resolution: float = SOURCE_MONITOR_LAMBDA_RESOLUTION,
    time_resolution: float = SOURCE_MONITOR_TIME_RESOLUTION,
) -> Monitor_nD:
    """Create the time and wavelength monitor at the entrance of W03-01-011.

    It is the entrance of W03-01-011: its opening, at the place and with the rotation of that
    guide on ``arm_iscs``, so it counts what enters the guide system. Its axes are those of
    :func:`.time_lambda_options` at ``GUIDE_START_Z``, without limits: the widest band the
    source simulates, ``L_MIN`` to ``L_MAX``, and the time window in which it reaches the
    monitor.

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
        A new monitor, relative to ``arm_iscs``.
    """
    to_gap_bbg, *_ = guides_bender()
    entrance = to_gap_bbg["guide_01_011"]
    return Monitor_nD(
        xwidth=GUIDE_01_011_IN_WIDTH / MM_PER_M,
        yheight=GUIDE_01_011_IN_HEIGHT / MM_PER_M,
        restore_neutron=1,
        options=time_lambda_options(GUIDE_START_Z, None, None, lambda_resolution, time_resolution),
        filename="source_monitor",  # the name of the monitor, in both models
        at=list(entrance.at),
        rotated=list(entrance.rotated),
        relative=entrance.relative,
    )


def _gap_monitor(
    name: str,
    *,
    width: float,
    height: float,
    path: float,
    z: float,
    relative: str,
    time_limits: tuple[float, float] | None,
    lambda_limits: tuple[float, float] | None,
    lambda_resolution: float,
    time_resolution: float,
) -> Monitor_nD:
    """Create a time and wavelength monitor in a gap of the guide.

    Parameters
    ----------
    name:
        Name of its output file, as the attribute of the instrument that holds it.
    width, height:
        Size of the opening, in mm.
    path:
        Flight path from the moderator, in m.
    z:
        Distance from the start of the gap, along it, in mm.
    relative:
        The arm at the start of the gap.
    time_limits, lambda_limits, lambda_resolution, time_resolution:
        As in :func:`.time_lambda_options`, at ``path``.

    Returns
    -------
    :
        A new monitor, relative to ``relative``.
    """
    return Monitor_nD(
        xwidth=width / MM_PER_M,
        yheight=height / MM_PER_M,
        restore_neutron=1,
        options=time_lambda_options(
            path, time_limits, lambda_limits, lambda_resolution, time_resolution
        ),
        filename=name,
        at=[0, 0, z / MM_PER_M],
        relative=relative,
    )


def beam_monitor(
    time_limits: tuple[float, float] | None = None,
    lambda_limits: tuple[float, float] | None = None,
    lambda_resolution: float = GAP_MONITOR_LAMBDA_RESOLUTION,
    time_resolution: float = GAP_MONITOR_TIME_RESOLUTION,
) -> Monitor_nD:
    """Create the time and wavelength monitor in the middle of the gap after the BBG, row 22.

    Parameters
    ----------
    time_limits, lambda_limits, lambda_resolution, time_resolution:
        As in :func:`.time_lambda_options`, at ``BEAM_MONITOR_PATH``.

    Returns
    -------
    :
        A new monitor, relative to ``arm_gap_bbg``.
    """
    return _gap_monitor(
        "beam_monitor",
        width=BEAM_MONITOR_WIDTH,
        height=BEAM_MONITOR_HEIGHT,
        path=BEAM_MONITOR_PATH,
        z=BEAM_MONITOR_Z,
        relative="arm_gap_bbg",
        time_limits=time_limits,
        lambda_limits=lambda_limits,
        lambda_resolution=lambda_resolution,
        time_resolution=time_resolution,
    )


def bwfa_monitor(
    time_limits: tuple[float, float] | None = None,
    lambda_limits: tuple[float, float] | None = None,
    lambda_resolution: float = GAP_MONITOR_LAMBDA_RESOLUTION,
    time_resolution: float = GAP_MONITOR_TIME_RESOLUTION,
) -> Monitor_nD:
    """Create the time and wavelength monitor in the middle of the gap after the BWFA, row 50.

    Parameters
    ----------
    time_limits, lambda_limits, lambda_resolution, time_resolution:
        As in :func:`.time_lambda_options`, at ``BWFA_MONITOR_PATH``.

    Returns
    -------
    :
        A new monitor, relative to ``arm_gap_bwfa``.
    """
    return _gap_monitor(
        "bwfa_monitor",
        width=BWFA_MONITOR_WIDTH,
        height=BWFA_MONITOR_HEIGHT,
        path=BWFA_MONITOR_PATH,
        z=BWFA_MONITOR_Z,
        relative="arm_gap_bwfa",
        time_limits=time_limits,
        lambda_limits=lambda_limits,
        lambda_resolution=lambda_resolution,
        time_resolution=time_resolution,
    )


def psc_monitor(
    time_limits: tuple[float, float] | None = None,
    lambda_limits: tuple[float, float] | None = None,
    lambda_resolution: float = GAP_MONITOR_LAMBDA_RESOLUTION,
    time_resolution: float = GAP_MONITOR_TIME_RESOLUTION,
) -> Monitor_nD:
    """Create the time and wavelength monitor after the pulse-shaping choppers.

    It is in their gap, halfway from ``ps2`` to the end of the gap.

    Parameters
    ----------
    time_limits, lambda_limits, lambda_resolution, time_resolution:
        As in :func:`.time_lambda_options`, at ``PSC_MONITOR_PATH``.

    Returns
    -------
    :
        A new monitor, relative to ``arm_psc``.
    """
    return _gap_monitor(
        "psc_monitor",
        width=PSC_MONITOR_WIDTH,
        height=PSC_MONITOR_HEIGHT,
        path=PSC_MONITOR_PATH,
        z=PSC_MONITOR_Z,
        relative="arm_psc",
        time_limits=time_limits,
        lambda_limits=lambda_limits,
        lambda_resolution=lambda_resolution,
        time_resolution=time_resolution,
    )


def mc_monitor(
    time_limits: tuple[float, float] | None = None,
    lambda_limits: tuple[float, float] | None = None,
    lambda_resolution: float = GAP_MONITOR_LAMBDA_RESOLUTION,
    time_resolution: float = GAP_MONITOR_TIME_RESOLUTION,
) -> Monitor_nD:
    """Create the time and wavelength monitor after the monochromatic choppers.

    It is in their gap, halfway from ``mc2`` to the end of the gap.

    Parameters
    ----------
    time_limits, lambda_limits, lambda_resolution, time_resolution:
        As in :func:`.time_lambda_options`, at ``MC_MONITOR_PATH``.

    Returns
    -------
    :
        A new monitor, relative to ``arm_mc``.
    """
    return _gap_monitor(
        "mc_monitor",
        width=MC_MONITOR_WIDTH,
        height=MC_MONITOR_HEIGHT,
        path=MC_MONITOR_PATH,
        z=MC_MONITOR_Z,
        relative="arm_mc",
        time_limits=time_limits,
        lambda_limits=lambda_limits,
        lambda_resolution=lambda_resolution,
        time_resolution=time_resolution,
    )


def sample_monitor(
    time_limits: tuple[float, float] | None = None,
    lambda_limits: tuple[float, float] | None = None,
    lambda_resolution: float = SAMPLE_MONITOR_LAMBDA_RESOLUTION,
    time_resolution: float = SAMPLE_MONITOR_TIME_RESOLUTION,
) -> Monitor_nD:
    """Create the time and wavelength monitor at the sample position.

    It records the beam reaching the sample. It is at ``sample_pos``, before the sample.

    Parameters
    ----------
    time_limits, lambda_limits, lambda_resolution, time_resolution:
        As in :func:`.time_lambda_options`, at ``SAMPLE_PATH``.

    Returns
    -------
    :
        A new monitor.
    """
    return Monitor_nD(
        xwidth=SAMPLE_MONITOR_WIDTH,
        yheight=SAMPLE_MONITOR_HEIGHT,
        restore_neutron=1,
        options=time_lambda_options(
            SAMPLE_PATH, time_limits, lambda_limits, lambda_resolution, time_resolution
        ),
        filename="sample_monitor",
        relative="sample_pos",
    )
