"""The choppers of T-Rex and their timing.

The disks take their frequency ``nu`` and their delay from the instrument parameters
``nu_<name>`` and ``delay_<name>``, which :func:`chopper_timings` computes from the central
wavelength, the repetition rate multiplication and the mode, as ``trex_primary.txt`` does.

The mode chooses the slit of the P- and M-Chopper disks that is on the beam: the wide one in
``"FF"`` (high flux), the narrow one in ``"RR"`` (high resolution). ``trex_primary.txt`` shifts the
slit edges for its high-flux mode; here the disks stay as they are and turn later, by the angle
from the slit to the beam, which is the same. The narrow slits are centered on ``beam_angle``,
so ``"RR"`` turns them by 0.

The flight path to a disk is its distance from the moderator along the beam axis: the position
of its gap in the Table of Optics plus the disk's place in the gap. The path through the bender
is taken to be that through guide 4a.
"""

from typing import Literal, NamedTuple, get_args

from mcdgs.core.components import (
    NXdisk_chopper,
)
from mcdgs.core.utils import speed_from_wavelength

from .guides.curve import ARM_BWC1_Z, ARM_BWC2_Z, BWC1_GAP_LENGTH, BWC2_GAP_LENGTH
from .guides.straight import ARM_MC_Z, ARM_PSC_Z, MC_GAP_LENGTH, PSC_GAP_LENGTH

Mode = Literal["FF", "RR"]  # high flux, high resolution


PSC_DISK_SPACING = 0.1  # FIXME, not in Table of Optics v3.2, in m
MC_DISK_SPACING = 0.01  # FIXME, not in Table of Optics v3.2, in m

# deg, opening and closing edge of each slit from the disk's zero mark
SLIT_EDGES = {
    "bw1": [-30.7, 30.7],
    "bw2": [-31.65, 31.65],
    "ps1": [-10, 10, 107.5, 142.5, 170, 190, 287.5, 322.5],
    "ps2": [-10, 10, 107.5, 142.5, 170, 190, 287.5, 322.5],
    "mc1": [178.75, 181.25, 362.8, 367.2],
    "mc2": [-1.25, 1.25, 182.8, 187.2],
}

# deg, where the beam crosses each disk from its zero mark, zero means disk axel below guide, 180 means above
BEAM_ANGLE = {"bw1": 0, "bw2": 0, "ps1": 180, "ps2": 180, "mc1": 180, "mc2": 0}

# index of the slit of each disk on the beam in each mode: wide for flux, narrow for resolution
MODE_SLIT: dict[Mode, dict[str, int]] = {
    "FF": {"bw1": 0, "bw2": 0, "ps1": 1, "ps2": 1, "mc1": 1, "mc2": 1},
    "RR": {"bw1": 0, "bw2": 0, "ps1": 2, "ps2": 2, "mc1": 0, "mc2": 0},
}

# m, distance of each disk from the entrance of its gap
DISK_Z = {
    "bw1": BWC1_GAP_LENGTH / 2,  # 0.015
    "bw2": BWC2_GAP_LENGTH / 2,  # 0.015
    "ps1": (PSC_GAP_LENGTH - PSC_DISK_SPACING) / 2,  # 0.05
    "ps2": (PSC_GAP_LENGTH + PSC_DISK_SPACING) / 2,  # 0.15
    "mc1": (MC_GAP_LENGTH - MC_DISK_SPACING) / 2,  # 0.0425
    "mc2": (MC_GAP_LENGTH + MC_DISK_SPACING) / 2,  # 0.0525
}

# m, distance of each disk from the moderator along the beam axis: its flight path
DISK_DISTANCE = {
    "bw1": ARM_BWC1_Z + DISK_Z["bw1"],  # 31.96404
    "bw2": ARM_BWC2_Z + DISK_Z["bw2"],  # 39.95504
    "ps1": ARM_PSC_Z + DISK_Z["ps1"],  # 107.95
    "ps2": ARM_PSC_Z + DISK_Z["ps2"],  # 108.05
    "mc1": ARM_MC_Z + DISK_Z["mc1"],  # 162.0275
    "mc2": ARM_MC_Z + DISK_Z["mc2"],  # 162.0375
}


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


def slit_angle(name: str, mode: Mode) -> float:
    """Compute the angle the disk turns past its zero mark to put the slit of a mode on the beam.

    Parameters
    ----------
    name:
        Chopper name.
    mode:
        ``"FF"`` or ``"RR"``.

    Returns
    -------
    :
        Angle from the center of the slit to ``beam_angle``, in [0, 360) deg.

    Examples
    --------
    >>> slit_angle("ps1", "FF"), slit_angle("mc1", "FF"), slit_angle("mc1", "RR")
    (55.0, 175.0, 0.0)
    """
    edges = SLIT_EDGES[name]
    slit = MODE_SLIT[mode][name]
    center = (edges[2 * slit] + edges[2 * slit + 1]) / 2
    return (BEAM_ANGLE[name] - center) % 360


def chopper_timings(
    lambda0: float, rrm: int, source_frequency: float, mode: Mode, t0_offset: float
) -> dict[str, Timing]:
    """Compute the timing of every chopper.

    The bandwidth choppers turn at the source frequency, the M-Chopper at ``rrm`` times it and
    the P-Chopper at ``PSC_FREQUENCY_RATIO`` of that, its second disk counter-rotating. Each disk
    has the slit of ``mode`` on the beam when a neutron of wavelength ``lambda0``, leaving at
    ``t0_offset``, reaches it.

    PSC_FREQUENCY_RATIO is 3/4 because the ratio of the distances from the source to the M-Chopper and from source to the P-Chopper is 162/108=3/2, and the P-Choppers has two openings per rotation.

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
    PSC_FREQUENCY_RATIO = 0.75

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
    for name, distance in DISK_DISTANCE.items():
        arrival = t0_offset + distance / speed
        # the zero mark passes the beam earlier, for the disk to have turned to the slit
        timings[name] = Timing(nu[name], arrival - slit_angle(name, mode) / (360 * nu[name]))
    return timings


def choppers() -> dict[str, NXdisk_chopper]:
    """Create the choppers.

    Returns
    -------
    :
        Choppers by name, in beam order, each relative to the arm of its gap.
    """
    return {
        "bw1": NXdisk_chopper(
            slit_edges=SLIT_EDGES["bw1"],
            radius=0.35,
            yheight=0.085,  # FIXME slit height unknown
            xwidth=0.0,  # 0.06
            nu="nu_bw1",
            delay="delay_bw1",
            beam_angle=BEAM_ANGLE["bw1"],
            zero_angle=0,
            abs_out=0,
            verbose=1,
            at=[0, 0, DISK_Z["bw1"]],
            relative="arm_bwc1",
        ),
        "bw2": NXdisk_chopper(
            slit_edges=SLIT_EDGES["bw2"],
            radius=0.35,
            yheight=0.085,  # FIXME slit height unknown
            xwidth=0.0,  # 0.06
            nu="nu_bw2",
            delay="delay_bw2",
            beam_angle=BEAM_ANGLE["bw2"],
            zero_angle=0,
            abs_out=0,
            verbose=1,
            at=[0, 0, DISK_Z["bw2"]],
            relative="arm_bwc2",
        ),
        "ps1": NXdisk_chopper(
            slit_edges=SLIT_EDGES["ps1"],
            radius=0.35,
            yheight=0.095,
            xwidth=0.0,  # 0.06
            nu="nu_ps1",
            delay="delay_ps1",
            beam_angle=BEAM_ANGLE["ps1"],
            zero_angle=0,
            abs_out=0,
            verbose=1,
            at=[0, 0, DISK_Z["ps1"]],
            relative="arm_psc",
        ),
        "ps2": NXdisk_chopper(
            slit_edges=SLIT_EDGES["ps2"],
            radius=0.35,
            yheight=0.095,
            xwidth=0.0,  # 0.06
            nu="nu_ps2",
            delay="delay_ps2",
            beam_angle=BEAM_ANGLE["ps2"],
            zero_angle=0,
            abs_out=0,
            verbose=1,
            at=[0, 0, DISK_Z["ps2"]],
            relative="arm_psc",
        ),
        "mc1": NXdisk_chopper(
            slit_edges=SLIT_EDGES["mc1"],
            radius=0.35,
            yheight=0.045,
            xwidth=0.0,  # 0.06
            nu="nu_mc1",
            delay="delay_mc1",
            beam_angle=BEAM_ANGLE["mc1"],
            zero_angle=0,
            abs_out=0,
            verbose=1,
            at=[0, 0, DISK_Z["mc1"]],
            relative="arm_mc",
        ),
        "mc2": NXdisk_chopper(
            slit_edges=SLIT_EDGES["mc2"],
            radius=0.35,
            yheight=0.045,
            xwidth=0.0,  # 0.06
            nu="nu_mc2",
            delay="delay_mc2",
            beam_angle=BEAM_ANGLE["mc2"],
            zero_angle=0,
            abs_out=0,
            verbose=1,
            at=[0, 0, DISK_Z["mc2"]],
            relative="arm_mc",
        ),
    }
