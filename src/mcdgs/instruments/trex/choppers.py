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

:func:`chopper_train` describes the six disks to chopper-lib, for the masked source.
"""

import math
from typing import Literal, NamedTuple, get_args

from mcdgs.core.components import Guide_gravity, NXdisk_chopper
from mcdgs.core.components.choppers import ChopperTrainParameters
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

# deg, opening and closing edge of each slit from the disk's zero mark
SLIT_EDGES = {
    "bw1": [-30.7, 30.7],
    "bw2": [-31.65, 31.65],
    "ps1": [-10, 10, 107.5, 142.5, 170, 190, 287.5, 322.5],
    "ps2": [-10, 10, 107.5, 142.5, 170, 190, 287.5, 322.5],
    "mc1": [178.75, 181.25, 362.8, 367.2],
    "mc2": [-1.25, 1.25, 182.8, 187.2],
}

# deg, where the beam crosses each disk from its zero mark: 0 puts the axle below the guide,
# 180 above it
BEAM_ANGLE = {"bw1": 0, "bw2": 0, "ps1": 180, "ps2": 180, "mc1": 180, "mc2": 0}

# m, radial height of the slits of each disk
SLIT_HEIGHT = {
    "bw1": 0.085,  # FIXME slit height unknown
    "bw2": 0.085,  # FIXME slit height unknown
    "ps1": 0.095,
    "ps2": 0.095,
    "mc1": 0.045,
    "mc2": 0.045,
}

# index of the slit of each disk on the beam in each mode: wide for flux, narrow for resolution
MODE_SLIT: dict[Mode, dict[str, int]] = {
    "FF": {"bw1": 0, "bw2": 0, "ps1": 1, "ps2": 1, "mc1": 1, "mc2": 1},
    "RR": {"bw1": 0, "bw2": 0, "ps1": 2, "ps2": 2, "mc1": 0, "mc2": 0},
}

# the arm at the entrance of the gap of each disk, which places it
GAP_ARM = {
    "bw1": "arm_bwc1",
    "bw2": "arm_bwc2",
    "ps1": "arm_psc",
    "ps2": "arm_psc",
    "mc1": "arm_mc",
    "mc2": "arm_mc",
}

# the guide in front of each disk, on the source side: the beam window is as wide as its exit
UPSTREAM_GUIDE = {
    "bw1": "guide_13a",
    "bw2": "guide_15a",
    "ps1": "guide_32",
    "ps2": "guide_32",
    "mc1": "guide_49",
    "mc2": "guide_49",
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

    The bandwidth choppers turn at the source frequency, the monochromatic choppers at ``rrm``
    times it and the pulse-shaping choppers at ``PSC_FREQUENCY_RATIO`` of that. Both pairs are
    counter-rotating: ``ps2`` turns the other way, and ``mc1`` and ``mc2`` turn the same way on
    opposite sides of the beam. Each disk
    has the slit of ``mode`` on the beam when a neutron of wavelength ``lambda0``, leaving at
    ``t0_offset``, reaches it.

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


def window_width(name: str) -> float:
    """Return the width of the beam window of a disk: the exit width of the guide in front.

    Parameters
    ----------
    name:
        Chopper name, such as ``"ps1"``: a key of :data:`UPSTREAM_GUIDE`.

    Returns
    -------
    :
        ``w2`` of the guide :data:`UPSTREAM_GUIDE` names, in m.

    Examples
    --------
    >>> round(window_width("mc1"), 4)
    0.0224
    """
    guides = {
        guide_name: guide
        for part in (*guides_curve(), *guides_straight())
        for guide_name, guide in part.items()
        if isinstance(guide, Guide_gravity)
    }
    return float(guides[UPSTREAM_GUIDE[name]].w2)


def chopper(name: str) -> NXdisk_chopper:
    """Create a chopper.

    Parameters
    ----------
    name:
        Chopper name, such as ``"ps1"``: a key of :data:`SLIT_EDGES`.

    Returns
    -------
    :
        A new chopper, relative to the arm of its gap, with its frequency and delay given by
        the instrument parameters ``nu_<name>`` and ``delay_<name>``, and the beam window of
        :func:`window_width`.
    """
    return NXdisk_chopper(
        slit_edges=SLIT_EDGES[name],
        radius=0.35,
        yheight=SLIT_HEIGHT[name],
        xwidth=window_width(name),
        nu=f"nu_{name}",
        delay=f"delay_{name}",
        beam_angle=BEAM_ANGLE[name],
        zero_angle=0,
        abs_out=0,
        verbose=1,
        at=[0, 0, DISK_Z[name]],
        relative=GAP_ARM[name],
    )


def aperture(disk: NXdisk_chopper) -> float:
    """Return the angular width of the beam on a disk, about its spindle, as chopper-lib takes it.

    The inner corners of the beam window, nearest the spindle, are the furthest round the disk.
    NXdisk_chopper hangs the window from the height of the rim at its edges,
    ``sqrt(radius**2 - (xwidth/2)**2)``, so those corners are ``yheight`` inside it.

    Parameters
    ----------
    disk:
        The disk, with its beam window.

    Returns
    -------
    :
        ``2 atan2(xwidth / 2, sqrt(radius**2 - (xwidth / 2)**2) - yheight)``, in deg.

    Examples
    --------
    >>> round(aperture(chopper("mc1")), 4)
    4.2148
    """
    half_width = disk.xwidth / 2
    reach = math.sqrt(disk.radius**2 - half_width**2) - disk.yheight
    return math.degrees(2 * math.atan2(half_width, reach))


def chopper_train() -> list[ChopperTrainParameters]:
    """Describe the six disks to chopper-lib, as the masked source takes them.

    Each disk is the one :func:`chopper` creates: the same edges, beam angle and window, and the
    same instrument parameters ``nu_<name>`` and ``delay_<name>``. Its path is
    :data:`DISK_DISTANCE`. Both McStas components put the disk angle ``a`` on the beam axis at
    ``delay + (beam_angle - a) / (360 nu)``, so the delay is the same.

    Returns
    -------
    :
        One entry per disk, in beam order.
    """
    train: list[ChopperTrainParameters] = []
    for name, path in DISK_DISTANCE.items():
        disk = chopper(name)
        train.append(
            {
                "speed": disk.nu,
                "delay": disk.delay,
                "beam": BEAM_ANGLE[name],
                "edges": list(disk.slit_edges),
                "path": path,
                "aperture": aperture(disk),
            }
        )
    return train
