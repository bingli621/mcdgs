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
side, which the instrument gives :func:`chopper`: the beam that reaches the disk. NXdisk_chopper hangs the window from the height of the
rim at its edges, so the slit height holds across the whole beam, and no neutron passes the
corners of the window outside the disk unchopped.

The flight path to a disk is its distance from the moderator along the beam axis: the position
of its gap in the Table of Optics plus the disk's place in the gap. The path through the bender
is taken to be that through guide 4a.

:data:`DISKS` describes the design of each disk as a :class:`TRexDisk`, with its flight path
and the slit of each mode.
"""

from dataclasses import dataclass
from typing import Literal, get_args

from mcdgs.core.components import Guide_gravity, NXdisk_chopper
from mcdgs.core.components.choppers import Disk
from mcdgs.core.utils import speed_from_wavelength

from .guides.curve import (
    ARM_BWC1_Z,
    ARM_BWC2_Z,
    BWC1_GAP_LENGTH,
    BWC2_GAP_LENGTH,
)
from .guides.straight import (
    ARM_MC_Z,
    ARM_PSC_Z,
    MC_GAP_LENGTH,
    PSC_GAP_LENGTH,
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

    def phase_from_mode(self, mode: Mode) -> float:
        """Compute the angle the disk turns past its zero mark to put the slit of a mode on the beam.

        Parameters
        ----------
        mode:
            ``"FF"`` or ``"RR"``.

        Returns
        -------
        :
            Angle from the center of the slit ``mode_slit`` gives for ``mode`` to
            ``beam_angle``, in [0, 360) deg.

        Examples
        --------
        >>> DISKS["ps1"].phase_from_mode("FF"), DISKS["mc1"].phase_from_mode("FF")
        (55.0, 175.0)
        >>> DISKS["mc1"].phase_from_mode("RR"), DISKS["mc2"].phase_from_mode("FF")
        (0.0, 185.0)
        """
        center, _ = self.slits[self.mode_slit[mode]]
        return float((self.beam_angle - center) % 360)


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
    ),
}


def chopper_timings(
    lambda0: float, rrm: int, source_frequency: float, mode: Mode, t0_offset: float
) -> dict[str, tuple[float, float]]:
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
        Signed rotation frequency in Hz, and time at which the disk's zero mark is on the beam
        in s, by chopper name, in beam order.

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
        phase = disk.phase_from_mode(mode)
        timings[name] = (nu[name], arrival - phase / (360 * nu[name]))
    return timings


def chopper(name: str, upstream_guide: Guide_gravity | None) -> NXdisk_chopper:
    """Create a chopper.

    Parameters
    ----------
    name:
        Chopper name, such as ``"ps1"``: a key of :data:`DISKS`.
    upstream_guide:
        The guide in front of the disk, on the source side: the beam window is as wide as its
        exit. ``None`` leaves the window unrestricted.

    Returns
    -------
    :
        A new chopper, relative to the arm of its gap, with its frequency and delay given by
        the instrument parameters ``nu_<name>`` and ``delay_<name>``, and a beam window of
        ``upstream_guide.w2``, or 0 without a guide.
    """
    disk = DISKS[name]
    return NXdisk_chopper(
        slit_edges=list(disk.slit_edges),
        radius=disk.radius,
        yheight=disk.slit_height,
        xwidth=0.0 if upstream_guide is None else float(upstream_guide.w2),
        nu=f"nu_{name}",
        delay=f"delay_{name}",
        beam_angle=disk.beam_angle,
        zero_angle=disk.zero_angle,
        abs_out=0,
        verbose=1,
        at=[0, 0, disk.z],
        relative=disk.gap_arm,
    )
