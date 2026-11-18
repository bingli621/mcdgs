"""The choppers of CSPEC.

:func:`chopper` creates one disk by name. The disks take their frequency ``nu`` and their delay
from the instrument parameters ``nu_<name>`` and ``delay_<name>``.

The specifications give where the first disk of each chopper turns: its distance from the start
of its gap, "gap before" in the Definitions sheet. The second disks of the pulse-shaping and the
monochromatic choppers, which make counter-rotating pairs with the first, are not in the sheet.
Nor are the disks themselves: each has one slit, centered on its zero mark, with the beam at the
top of the disk.

Each disk has a beam window as wide as the exit of the guide in front of it, on the source
side, which the instrument gives :func:`chopper`.

The flight path to a disk is the trajectory from the moderator to the start of its gap, from
the sheet, plus the disk's place in the gap.

:data:`DISKS` describes the design of each disk as a :class:`.Disk`, with its flight path.
"""

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
PSC_DISK_SPACING = 6  # FIXME, not in CSPEC nl-specs 1.E, in mm
MC_DISK_SPACING = 6  # FIXME, not in CSPEC nl-specs 1.E, in mm

DISK_RADIUS = 0.35  # FIXME, not in CSPEC nl-specs 1.E, in m


# the seven disks, in beam order
DISKS: dict[str, Disk] = {
    "bw1": Disk(
        slits=((0, 40.7),),  # FIXME, not in CSPEC nl-specs 1.E
        radius=DISK_RADIUS,
        beam_angle=0,
        zero_angle=0,
        slit_height=0.055,  # FIXME, not in CSPEC nl-specs 1.E
        gap_arm="arm_bwc1",
        gap_z=ARM_BWC1_PATH,
        z=BWC1_GAP_BEFORE / MM_PER_M,  # distance 14.93342
    ),
    "bw2": Disk(
        slits=((0, 41.9),),  # FIXME, not in CSPEC nl-specs 1.E
        radius=DISK_RADIUS,
        beam_angle=0,
        zero_angle=0,
        slit_height=0.055,  # FIXME, not in CSPEC nl-specs 1.E
        gap_arm="arm_bwc2",
        gap_z=ARM_BWC2_PATH,
        z=BWC2_GAP_BEFORE / MM_PER_M,  # distance 20.45545
    ),
    "bw3": Disk(
        slits=((0, 193),),  # FIXME, not in CSPEC nl-specs 1.E
        radius=DISK_RADIUS,
        beam_angle=0,
        zero_angle=0,
        slit_height=0.0739,  # FIXME, not in CSPEC nl-specs 1.E
        gap_arm="arm_bwc3",
        gap_z=ARM_BWC3_PATH,
        z=BWC3_GAP_BEFORE / MM_PER_M,  # distance 104.61213
    ),
    "ps1": Disk(
        slits=((0, 24.23),),  # FIXME, not in CSPEC nl-specs 1.E
        radius=DISK_RADIUS,
        beam_angle=0,
        zero_angle=0,
        slit_height=0.07389,  # FIXME, not in CSPEC nl-specs 1.E
        gap_arm="arm_psc",
        gap_z=ARM_PSC_PATH,
        z=PSC_GAP_BEFORE / MM_PER_M,  # distance 105.66872
    ),
    "ps2": Disk(
        slits=((0, 24.23),),  # FIXME, not in CSPEC nl-specs 1.E
        radius=DISK_RADIUS,
        beam_angle=0,
        zero_angle=0,
        slit_height=0.07389,  # FIXME, not in CSPEC nl-specs 1.E
        gap_arm="arm_psc",
        gap_z=ARM_PSC_PATH,
        z=(PSC_GAP_BEFORE + PSC_DISK_SPACING) / MM_PER_M,  # distance 105.67472
    ),
    "mc1": Disk(
        slits=((0, 4.45),),  # FIXME, not in CSPEC nl-specs 1.E
        radius=DISK_RADIUS,
        beam_angle=0,
        zero_angle=0,
        slit_height=0.049,  # FIXME, not in CSPEC nl-specs 1.E
        gap_arm="arm_mc",
        gap_z=ARM_MC_PATH,
        z=MC_GAP_BEFORE / MM_PER_M,  # distance 158.50404
    ),
    "mc2": Disk(
        slits=((0, 4.45),),  # FIXME, not in CSPEC nl-specs 1.E
        radius=DISK_RADIUS,
        beam_angle=0,
        zero_angle=0,
        slit_height=0.049,  # FIXME, not in CSPEC nl-specs 1.E
        gap_arm="arm_mc",
        gap_z=ARM_MC_PATH,
        z=(MC_GAP_BEFORE + MC_DISK_SPACING) / MM_PER_M,  # distance 158.51004
    ),
}


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
        A new chopper, relative to the arm of its gap, with its frequency and delay given by
        the instrument parameters ``nu_<name>`` and ``delay_<name>``, and a beam window of
        ``upstream_guide.w2``, or 0 without a guide; and the disk's flight path from the
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
        beam_angle=disk.beam_angle,
        zero_angle=disk.zero_angle,
        abs_out=0,
        verbose=1,
        at=[0, 0, disk.z],
        relative=disk.gap_arm,
    )
    return nx_chopper, disk.distance
