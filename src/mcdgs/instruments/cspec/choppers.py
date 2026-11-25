"""The choppers of CSPEC and their timing.

:func:`chopper` creates one disk by name. The disks take their frequency ``nu``, their delay and
their ``park_angle`` from the instrument parameters ``nu_<name>``, ``delay_<name>`` and
``park_<name>``, which :func:`chopper_timings` computes from the central wavelength and the
frequencies of the monochromatic and the RRM choppers. A frequency of 0 of the monochromatic
choppers parks every disk on its first slit.

CSPEC_nl-specs_1.E.xlsx gives where the first disk of the pulse-shaping and of the monochromatic
choppers turns, and where each bandwidth chopper turns: the distance from the start of its gap,
"gap before" in the Definitions sheet. The second disk of each counter-rotating pair, the RRM
chopper and the disks themselves are not in the sheet.

Each disk's slits are centered on angles from its zero mark. A disk whose axle is above the beam
has it at ``beam_angle`` 180, one whose axle is below at 0. :meth:`CSPECDisk.phase` is the
angle the disk turns past its zero mark to put the center of its first slit on the beam.

Each disk has a beam window as wide as the exit of the guide in front of it, on the source
side, which the instrument gives :func:`chopper`.

The flight path to a disk is the trajectory from the moderator to the start of its gap, from
the sheet, plus the disk's place in the gap.

:data:`DISKS` describes the design of each disk as a :class:`CSPECDisk`, with its flight path.
"""

from dataclasses import dataclass
from typing import NamedTuple

from mcdgs.core import speed_from_wavelength
from mcdgs.core.components import Guide_gravity, NXdisk_chopper
from mcdgs.core.components.choppers import Disk

from .guides.bender import ARM_BWC1_PATH, ARM_BWC2_PATH
from .guides.frame import MM_PER_M
from .guides.straight import ARM_BWC3_PATH, ARM_MC_PATH, ARM_PSC_PATH

BWC1_GAP_BEFORE = 15  # mm, cell E12 of Definitions
BWC2_GAP_BEFORE = 15  # mm, cell E16 of Definitions
BWC3_GAP_BEFORE = 15  # mm, cell E24 of Definitions
PSC_GAP_BEFORE = 12.71  # mm, cell E28 of Definitions
MC_GAP_BEFORE = 48.22  # mm, cell E32 of Definitions

# from the first disk of a pair to the second
PSC_DISK_SPACING = 10  # FIXME, not in CSPEC_nl-specs_1.E.xlsx, in mm
MC_DISK_SPACING = 10  # FIXME, not in CSPEC_nl-specs_1.E.xlsx, in mm
# the RRM chopper turns in the gap of the monochromatic choppers, halfway from its start to mc1.
# Not at the start itself: the sheet has it 27 nm before the exit of W03-09-01, so a disk there
# would be behind the neutrons leaving the guide
RRM_GAP_BEFORE = MC_GAP_BEFORE / 2  # FIXME, not in CSPEC_nl-specs_1.E.xlsx, in mm

DISK_RADIUS = 0.35  # FIXME, not in CSPEC_nl-specs_1.E.xlsx, in m
SLIT_HEIGHT = 0.1  # FIXME, not in CSPEC_nl-specs_1.E.xlsx, in m

# the frequency of the pulse-shaping choppers over that of the monochromatic ones, see
# chopper_timings. FIXME, not in CSPEC_nl-specs_1.E.xlsx
PSC_FREQUENCY_RATIO = 0.5


@dataclass(frozen=True, kw_only=True)
class CSPECDisk(Disk):
    """The design of one chopper disk of CSPEC: a :class:`.Disk` phased on its first slit."""

    def phase(self) -> float:
        """Compute the angle the disk turns past its zero mark to put its first slit on the beam.

        Returns
        -------
        :
            Angle from the center of the first slit to ``beam_angle``, in [0, 360) deg.

        Examples
        --------
        >>> DISKS["bw1"].phase(), DISKS["mc2"].phase()
        (180.0, 0.0)
        """
        center, _ = self.slits[0]
        return float((self.beam_angle - center) % 360)


# the eight disks, in beam order
DISKS: dict[str, CSPECDisk] = {
    "bw1": CSPECDisk(
        slits=((0, 40.7),),  # FIXME, not in CSPEC_nl-specs_1.E.xlsx
        radius=DISK_RADIUS,
        beam_angle=180,  # FIXME, not in CSPEC_nl-specs_1.E.xlsx
        zero_angle=0,
        slit_height=SLIT_HEIGHT,
        gap_arm="arm_bwc1",
        gap_z=ARM_BWC1_PATH,
        z=BWC1_GAP_BEFORE / MM_PER_M,  # distance 14.93342
    ),
    "bw2": CSPECDisk(
        slits=((0, 41.7),),  # FIXME, not in CSPEC_nl-specs_1.E.xlsx
        radius=DISK_RADIUS,
        beam_angle=180,  # FIXME, not in CSPEC_nl-specs_1.E.xlsx
        zero_angle=0,
        slit_height=SLIT_HEIGHT,
        gap_arm="arm_bwc2",
        gap_z=ARM_BWC2_PATH,
        z=BWC2_GAP_BEFORE / MM_PER_M,  # distance 20.45545
    ),
    "bw3": CSPECDisk(
        slits=((0, 193),),  # FIXME, not in CSPEC_nl-specs_1.E.xlsx
        radius=DISK_RADIUS,
        beam_angle=180,  # FIXME, not in CSPEC_nl-specs_1.E.xlsx
        zero_angle=0,
        slit_height=SLIT_HEIGHT,
        gap_arm="arm_bwc3",
        gap_z=ARM_BWC3_PATH,
        z=BWC3_GAP_BEFORE / MM_PER_M,  # distance 104.61213
    ),
    "ps1": CSPECDisk(
        slits=((0, 24.23), (120, 24.23), (240, 24.23)),  # FIXME, not in CSPEC_nl-specs_1.E.xlsx
        radius=DISK_RADIUS,
        beam_angle=180,  # FIXME, not in CSPEC_nl-specs_1.E.xlsx
        zero_angle=0,
        slit_height=SLIT_HEIGHT,
        gap_arm="arm_psc",
        gap_z=ARM_PSC_PATH,
        z=PSC_GAP_BEFORE / MM_PER_M,  # distance 105.66872
    ),
    "ps2": CSPECDisk(
        slits=((0, 24.23), (120, 24.23), (240, 24.23)),  # FIXME, not in CSPEC_nl-specs_1.E.xlsx
        radius=DISK_RADIUS,
        beam_angle=0,  # FIXME, not in CSPEC_nl-specs_1.E.xlsx
        zero_angle=0,
        slit_height=SLIT_HEIGHT,
        gap_arm="arm_psc",
        gap_z=ARM_PSC_PATH,
        z=(PSC_GAP_BEFORE + PSC_DISK_SPACING) / MM_PER_M,  # distance 105.67872
    ),
    "rrm": CSPECDisk(
        slits=((0, 4.45),),  # FIXME, not in CSPEC_nl-specs_1.E.xlsx
        radius=DISK_RADIUS,
        beam_angle=0,  # FIXME, not in CSPEC_nl-specs_1.E.xlsx
        zero_angle=0,
        slit_height=SLIT_HEIGHT,
        gap_arm="arm_mc",
        gap_z=ARM_MC_PATH,
        z=RRM_GAP_BEFORE / MM_PER_M,  # distance 158.47993
    ),
    "mc1": CSPECDisk(
        slits=((0, 4.45),),  # FIXME, not in CSPEC_nl-specs_1.E.xlsx
        radius=DISK_RADIUS,
        beam_angle=180,  # FIXME, not in CSPEC_nl-specs_1.E.xlsx
        zero_angle=0,
        slit_height=SLIT_HEIGHT,
        gap_arm="arm_mc",
        gap_z=ARM_MC_PATH,
        z=MC_GAP_BEFORE / MM_PER_M,  # distance 158.50404
    ),
    "mc2": CSPECDisk(
        slits=((0, 4.45),),  # FIXME, not in CSPEC_nl-specs_1.E.xlsx
        radius=DISK_RADIUS,
        beam_angle=0,  # FIXME, not in CSPEC_nl-specs_1.E.xlsx
        zero_angle=0,
        slit_height=SLIT_HEIGHT,
        gap_arm="arm_mc",
        gap_z=ARM_MC_PATH,
        z=(MC_GAP_BEFORE + MC_DISK_SPACING) / MM_PER_M,  # distance 158.51404
    ),
}


class ChopperTiming(NamedTuple):
    """The timing of one chopper, as the instrument parameters of its ``NXdisk_chopper``.

    Attributes
    ----------
    nu:
        Signed rotation frequency, in Hz. 0 parks the disk.
    delay:
        Time at which the disk's zero mark is on the beam, in s. 0 while parked.
    park_angle:
        Where the zero mark stands while parked, in deg: on the first slit. 0 while turning,
        where ``NXdisk_chopper`` does not read it.
    """

    nu: float
    delay: float
    park_angle: float


def chopper_timings(
    lambda0: float | None,
    nu_mc: float,
    n: int,
    source_frequency: float,
    t0_offset: float,
    bw3_delay: float = 0.0,
) -> dict[str, ChopperTiming]:
    """Compute the timing of every chopper.

    The bandwidth choppers turn at the source frequency, the monochromatic choppers at ``nu_mc``,
    the pulse-shaping choppers at ``PSC_FREQUENCY_RATIO`` of it and the RRM chopper at
    ``nu_mc / n``. Both pairs are counter-rotating, their slit edges crossing the beam in opposite
    directions: ``ps2`` and ``mc2`` turn the same way as ``ps1`` and ``mc1``, from the other side of
    the beam, their axles below it where those of ``ps1`` and ``mc1`` are above. Each
    disk has its first slit centered on the beam when a neutron of wavelength ``lambda0``,
    leaving at ``t0_offset``, reaches it. ``bw3`` turns ``bw3_delay`` later still.

    ``nu_mc`` 0 parks every disk, the bandwidth choppers too: each stands still with its first
    slit on the beam, its zero mark turned by :meth:`CSPECDisk.phase`. Nothing is then phased to
    a wavelength, so ``lambda0`` may be ``None``, and ``bw3_delay`` is not used.

    Parameters
    ----------
    lambda0:
        Central wavelength, in Å. ``None`` only with ``nu_mc`` 0.
    nu_mc:
        Frequency of the monochromatic choppers, in Hz. 0 parks every chopper.
    n:
        The RRM chopper turns at ``nu_mc / n``.
    source_frequency:
        Pulse frequency of the source, in Hz.
    t0_offset:
        Time from the start of the pulse at which the neutrons the choppers are phased to leave
        the moderator, in s.
    bw3_delay:
        Delay added to that of ``bw3``, in s.

    Returns
    -------
    :
        The timing of each chopper by name, in beam order.

    Raises
    ------
    ValueError
        If ``nu_mc`` is negative, or if ``lambda0`` is ``None`` while ``nu_mc`` is not 0.

    Examples
    --------
    >>> chopper_timings(None, 0.0, 1, 14, 0.0017)["bw1"]
    ChopperTiming(nu=0.0, delay=0.0, park_angle=180.0)
    """
    if nu_mc < 0:
        raise ValueError(f"nu_mc {nu_mc} is negative")
    if nu_mc == 0:
        return {
            name: ChopperTiming(nu=0.0, delay=0.0, park_angle=disk.phase())
            for name, disk in DISKS.items()
        }
    if lambda0 is None:
        raise ValueError("lambda0 is None: the choppers have no wavelength to be phased to")
    nu_psc = PSC_FREQUENCY_RATIO * nu_mc
    nu = {
        "bw1": source_frequency,
        "bw2": source_frequency,
        "bw3": source_frequency,
        "ps1": nu_psc,
        "ps2": nu_psc,  # the same way, from below the beam: counter-rotating with ps1
        "rrm": nu_mc / n,
        "mc1": nu_mc,
        "mc2": nu_mc,  # the same way, from below the beam: counter-rotating with mc1
    }
    speed = speed_from_wavelength(lambda0)
    timings = {}
    for name, disk in DISKS.items():
        arrival = t0_offset + disk.distance / speed
        # the zero mark passes the beam earlier, for the disk to have turned to the slit
        delay = arrival - disk.phase() / (360 * nu[name])
        timings[name] = ChopperTiming(nu=nu[name], delay=delay, park_angle=0.0)
    timings["bw3"] = timings["bw3"]._replace(delay=timings["bw3"].delay + bw3_delay)
    return timings


def chopper(name: str, upstream_guide: Guide_gravity | None) -> tuple[NXdisk_chopper, float]:
    """Create a chopper, with its flight path.

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
        A new chopper, relative to the arm of its gap, with its frequency, delay and park
        angle given by the instrument parameters ``nu_<name>``, ``delay_<name>`` and
        ``park_<name>``, and a beam window of ``upstream_guide.w2``, or 0 without a guide; and the disk's flight path from the
        moderator, in m, which the chopper does not know.
    """
    disk = DISKS[name]
    nx_chopper = NXdisk_chopper(
        slit_edges=list(disk.slit_edges),
        radius=disk.radius,
        yheight=disk.slit_height,
        xwidth=0.0 if upstream_guide is None else float(upstream_guide.w2),
        nu=f"nu_{name}",
        delay=f"delay_{name}",
        park_angle=f"park_{name}",
        beam_angle=disk.beam_angle,
        zero_angle=disk.zero_angle,
        abs_out=0,
        verbose=1,
        at=[0, 0, disk.z],
        relative=disk.gap_arm,
    )
    return nx_chopper, disk.distance
