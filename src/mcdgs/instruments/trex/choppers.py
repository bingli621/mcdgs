"""The choppers of T-Rex and their timing.

:func:`chopper` creates one disk by name. The disks take their frequency ``nu`` and their delay
from the instrument parameters ``nu_<name>`` and ``delay_<name>``, which :func:`chopper_timings`
computes from the central wavelength, the repetition rate multiplication and the mode.

The mode chooses the slit of the pulse-shaping and monochromatic chopper disks that is on the
beam: the wide one in ``"FF"`` (high flux), the narrow one in ``"RR"`` (high resolution). Rather
than shift the slit edges, the disks stay as they are and turn later, by the angle from the slit
to the beam, which is the same. The narrow slits are centered on ``beam_angle``, so ``"RR"``
turns them by 0.

Each disk has a beam window as wide as the exit of the guide in front of it, on the source
side: the beam that reaches the disk. NXdisk_chopper hangs the window from the height of the
rim at its edges, so the slit height holds across the whole beam, and no neutron passes the
corners of the window outside the disk unchopped.

The flight path to a disk is its distance from the moderator along the beam axis: the position
of its gap in the Table of Optics plus the disk's place in the gap. The path through the bender
is taken to be that through guide 4a.

:data:`DISKS` describes the design of each disk as a :class:`TRexDisk`, with its flight path
and the slit of each mode.
"""

from dataclasses import dataclass
from typing import Literal, NamedTuple, get_args

from mcdgs.core.components import Guide_gravity, NXdisk_chopper
from mcdgs.core.components.choppers import Disk
from mcdgs.core.utils import speed_from_wavelength

from .guides.curve import (
    ARM_BWC1_Z,
    ARM_BWC2_Z,
    BWC1_GAP_LENGTH,
    BWC2_GAP_LENGTH,
    guides_curve,
)
from .guides.straight import (
    ARM_MC_Z,
    ARM_PSC_Z,
    MC_GAP_LENGTH,
    PSC_GAP_LENGTH,
    guides_straight,
)

Mode = Literal["FF", "RR"]  # high flux, high resolution


PSC_DISK_SPACING = 0.1  # FIXME, not in Table of Optics v3.2, in m
MC_DISK_SPACING = 0.01  # FIXME, not in Table of Optics v3.2, in m

# the frequency of the pulse-shaping choppers over that of the monochromatic ones, see
# chopper_timings
PSC_FREQUENCY_RATIO = 0.75


@dataclass(frozen=True, kw_only=True)
class TRexDisk(Disk):
    """The design of one chopper disk of T-Rex: a :class:`.Disk` with a slit for each mode.

    Attributes
    ----------
    mode_slit:
        Index in ``slits`` of the slit on the beam in each mode: wide for flux, narrow for
        resolution.
    """

    mode_slit: dict[Mode, int]


# the six disks, in beam order
DISKS: dict[str, TRexDisk] = {
    "bw1": TRexDisk(
        slits=((0, 61.4),),
        radius=0.35,
        beam_angle=0,
        zero_angle=0,
        slit_height=0.085,  # FIXME slit height unknown
        mode_slit={"FF": 0, "RR": 0},
        gap_arm="arm_bwc1",
        gap_z=ARM_BWC1_Z,
        z=BWC1_GAP_LENGTH / 2,  # 0.015, distance 31.96404
        upstream_guide="guide_13a",
    ),
    "bw2": TRexDisk(
        slits=((0, 63.3),),
        radius=0.35,
        beam_angle=0,
        zero_angle=0,
        slit_height=0.085,  # FIXME slit height unknown
        mode_slit={"FF": 0, "RR": 0},
        gap_arm="arm_bwc2",
        gap_z=ARM_BWC2_Z,
        z=BWC2_GAP_LENGTH / 2,  # 0.015, distance 39.95504
        upstream_guide="guide_15a",
    ),
    "ps1": TRexDisk(
        slits=((0, 20), (125, 35), (180, 20), (305, 35)),
        radius=0.35,
        beam_angle=180,
        zero_angle=0,
        slit_height=0.095,
        mode_slit={"FF": 1, "RR": 2},
        gap_arm="arm_psc",
        gap_z=ARM_PSC_Z,
        z=(PSC_GAP_LENGTH - PSC_DISK_SPACING) / 2,  # 0.05, distance 107.95
        upstream_guide="guide_32",
    ),
    "ps2": TRexDisk(
        slits=((0, 20), (125, 35), (180, 20), (305, 35)),
        radius=0.35,
        beam_angle=180,
        zero_angle=0,
        slit_height=0.095,
        mode_slit={"FF": 1, "RR": 2},
        gap_arm="arm_psc",
        gap_z=ARM_PSC_Z,
        z=(PSC_GAP_LENGTH + PSC_DISK_SPACING) / 2,  # 0.15, distance 108.05
        upstream_guide="guide_32",
    ),
    "mc1": TRexDisk(
        slits=((180, 2.5), (365, 4.4)),
        radius=0.35,
        beam_angle=180,
        zero_angle=0,
        slit_height=0.045,
        mode_slit={"FF": 1, "RR": 0},
        gap_arm="arm_mc",
        gap_z=ARM_MC_Z,
        z=(MC_GAP_LENGTH - MC_DISK_SPACING) / 2,  # 0.0425, distance 162.0275
        upstream_guide="guide_49",
    ),
    # the wide slit is centered on 175, not on the 185 of mctrex: the narrow slit then comes by
    # the beam 175 deg of turn after the wide one, and that of mc1 185 deg after, so the two are
    # never open together
    "mc2": TRexDisk(
        slits=((0, 2.5), (175, 4.4)),
        radius=0.35,
        beam_angle=0,
        zero_angle=0,
        slit_height=0.045,
        mode_slit={"FF": 1, "RR": 0},
        gap_arm="arm_mc",
        gap_z=ARM_MC_Z,
        z=(MC_GAP_LENGTH + MC_DISK_SPACING) / 2,  # 0.0525, distance 162.0375
        upstream_guide="guide_49",
    ),
}


def phase_shift(disk: Disk, angle: float) -> float:
    """Compute the angle a disk turns past its zero mark to put an angle on the beam.

    Parameters
    ----------
    disk:
        The disk.
    angle:
        Angle on the disk from its zero mark, such as the center of a slit, in deg.

    Returns
    -------
    :
        Angle from ``angle`` to the disk's ``beam_angle``, in [0, 360) deg.

    Examples
    --------
    >>> phase_shift(DISKS["ps1"], 125), phase_shift(DISKS["mc1"], 365)
    (55.0, 175.0)
    """
    return float((disk.beam_angle - angle) % 360)


def trex_guides() -> dict[str, Guide_gravity]:
    """Return the guides of T-Rex by name, in beam order.

    Returns
    -------
    :
        The :class:`.Guide_gravity` components of :func:`.guides_curve` and
        :func:`.guides_straight`.
    """
    return {
        guide_name: guide
        for part in (*guides_curve(), *guides_straight())
        for guide_name, guide in part.items()
        if isinstance(guide, Guide_gravity)
    }


def mode_phase_shift(name: str, mode: Mode) -> float:
    """Compute the phase shift that puts the center of the slit of a mode on the beam.

    Parameters
    ----------
    name:
        Chopper name, such as ``"ps1"``: a key of :data:`DISKS`.
    mode:
        ``"FF"`` or ``"RR"``.

    Returns
    -------
    :
        :func:`phase_shift` of the center of the disk's ``mode_slit`` of ``mode``, in deg.

    Examples
    --------
    >>> mode_phase_shift("ps1", "FF"), mode_phase_shift("mc1", "FF"), mode_phase_shift("mc1", "RR")
    (55.0, 175.0, 0.0)
    """
    disk = DISKS[name]
    center, _ = disk.slits[disk.mode_slit[mode]]
    return phase_shift(disk, center)


class Timing(NamedTuple):
    """The timing of a chopper.

    Attributes
    ----------
    nu:
        Signed rotation frequency, in Hz.
    delay:
        Time at which the disk's zero mark is on the beam, in s.
    """

    nu: float
    delay: float


def chopper_timings(
    lambda0: float, rrm: int, source_frequency: float, mode: Mode, t0_offset: float
) -> dict[str, Timing]:
    """Compute the timing of every chopper.

    The bandwidth choppers turn at the source frequency, the monochromatic choppers at ``rrm``
    times it and the pulse-shaping choppers at ``PSC_FREQUENCY_RATIO`` of that. Both pairs are
    counter-rotating: ``ps2`` turns the other way, and ``mc1`` and ``mc2`` turn the same way on
    opposite sides of the beam. Each disk has the slit of ``mode`` on the beam when a neutron of
    wavelength ``lambda0``, leaving at ``t0_offset``, reaches it.

    PSC_FREQUENCY_RATIO is 3/4 because the ratio of the distances from the source to the
    monochromatic choppers and to the pulse-shaping choppers is 162/108 = 3/2, and the
    pulse-shaping choppers have two openings per rotation.

    Parameters
    ----------
    lambda0:
        Central wavelength, in Å.
    rrm:
        Repetition rate multiplication factor.
    source_frequency:
        Pulse frequency of the source, in Hz.
    mode:
        ``"FF"`` for high flux, ``"RR"`` for high resolution.
    t0_offset:
        Time from the start of the pulse at which the neutrons the choppers are phased to leave
        the moderator, in s.

    Returns
    -------
    :
        Timing by chopper name, in beam order.

    Raises
    ------
    ValueError
        If ``mode`` is not ``"FF"`` or ``"RR"``.
    """
    if mode not in get_args(Mode):
        raise ValueError(f"mode {mode!r} is not 'FF' or 'RR'")
    nu_mc = rrm * source_frequency
    nu_psc = PSC_FREQUENCY_RATIO * nu_mc
    nu = {
        "bw1": source_frequency,
        "bw2": source_frequency,
        "ps1": nu_psc,
        "ps2": -nu_psc,
        "mc1": nu_mc,
        "mc2": nu_mc,
    }
    speed = speed_from_wavelength(lambda0)
    timings = {}
    for name, disk in DISKS.items():
        arrival = t0_offset + disk.distance / speed
        # the zero mark passes the beam earlier, for the disk to have turned to the slit
        shift = mode_phase_shift(name, mode)
        timings[name] = Timing(nu[name], arrival - shift / (360 * nu[name]))
    return timings


def chopper(name: str) -> NXdisk_chopper:
    """Create a chopper.

    Parameters
    ----------
    name:
        Chopper name, such as ``"ps1"``: a key of :data:`DISKS`.

    Returns
    -------
    :
        A new chopper, relative to the arm of its gap, with its frequency and delay given by
        the instrument parameters ``nu_<name>`` and ``delay_<name>``, and the beam window of
        :meth:`.Disk.window_width` in the guides of :func:`trex_guides`.
    """
    disk = DISKS[name]
    return NXdisk_chopper(
        slit_edges=list(disk.slit_edges),
        radius=disk.radius,
        yheight=disk.slit_height,
        xwidth=disk.window_width(trex_guides()),
        nu=f"nu_{name}",
        delay=f"delay_{name}",
        beam_angle=disk.beam_angle,
        zero_angle=disk.zero_angle,
        abs_out=0,
        verbose=1,
        at=[0, 0, disk.z],
        relative=disk.gap_arm,
    )
