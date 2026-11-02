"""The choppers of T-Rex."""

from mcdgs.core.components import (
    NXdisk_chopper,
)

from .guides.curve import BWC1_GAP_LENGTH, BWC2_GAP_LENGTH
from .guides.straight import MC_GAP_LENGTH, PSC_GAP_LENGTH

bw1 = NXdisk_chopper(
    slit_edges=[-30.7, 30.7],
    radius=0.35,
    yheight=0.085,  # FIXME
    xwidth=0.0,  # 0.06
    beam_angle=0,
    zero_angle=0,
    abs_out=0,
    verbose=1,
    at=[0, 0, BWC1_GAP_LENGTH / 2],
    relative="arm_bwc1",
)

bw2 = NXdisk_chopper(
    slit_edges=[-31.65, 31.65],
    radius=0.35,
    yheight=0.085,  # FIXME
    xwidth=0.0,  # 0.06
    beam_angle=0,
    zero_angle=0,
    abs_out=0,
    verbose=1,
    at=[0, 0, BWC2_GAP_LENGTH / 2],
    relative="arm_bwc2",
)

ps1 = NXdisk_chopper(
    slit_edges=[-10, 10, 107.5, 142.5, 170, 190, 287.5, 322.5],
    radius=0.35,
    yheight=0.095,
    xwidth=0.0,  # 0.06
    beam_angle=180,
    zero_angle=0,
    abs_out=0,
    verbose=1,
    at=[0, 0, PSC_GAP_LENGTH / 2 - 0.05],
    relative="arm_psc",
)

ps2 = NXdisk_chopper(
    slit_edges=[-10, 10, 107.5, 142.5, 170, 190, 287.5, 322.5],
    radius=0.35,
    yheight=0.095,
    xwidth=0.0,  # 0.06
    beam_angle=180,
    zero_angle=0,
    abs_out=0,
    verbose=1,
    at=[0, 0, PSC_GAP_LENGTH / 2 + 0.05],
    relative="arm_psc",
)

mc1 = NXdisk_chopper(
    slit_edges=[178.75, 181.25, 362.8, 367.2],
    radius=0.35,
    yheight=0.045,
    xwidth=0.0,  # 0.06
    beam_angle=180,
    zero_angle=0,
    abs_out=0,
    verbose=1,
    at=[0, 0, MC_GAP_LENGTH / 2 - 0.005],
    relative="arm_mc",
)

mc2 = NXdisk_chopper(
    slit_edges=[-1.25, 1.25, 182.8, 187.2],
    radius=0.35,
    yheight=0.045,
    xwidth=0.0,  # 0.06
    beam_angle=0,
    zero_angle=0,
    abs_out=0,
    verbose=1,
    at=[0, 0, MC_GAP_LENGTH / 2 + 0.005],
    relative="arm_mc",
)
