"""The T-Rex instrument and its control parameters."""

from dataclasses import dataclass, replace

from mcdgs.core.base import Instrument
from mcdgs.core.components import (
    Arm,
    ESS_butterfly,
    Guide_gravity,
    Pol_bender,
    Progress_bar,
    Single_crystal,
    Slit,
)
from mcdgs.core.params import InstrParam, param

from .guides.curve import guides_curve
from .guides.extraction import BENDER_ROTATION, guides_extraction
from .guides.straight import guides_straight, sample_pos, slit_mc
from .sources import source


@dataclass(frozen=True)
class Params:
    """Control parameters of T-Rex.

    Frozen: change them with ``TRex.update_params``. Each field's unit and comment are in its
    metadata.

    Parameters
    ----------
    lambda0:
        Central wavelength, in Å.
    rrm:
        Repetition rate multiplication factor.
    omega:
        Goniometer rotation about +y, in deg.
    bender_in:
        True uses the bender and the cold guide 4b, False uses the thermal guide 4a.
    b_rot:
        Rotation of the bender about +y, in deg.

    """

    lambda0: float = param(2.0, "Å", "central wavelength")
    rrm: int = param(18, "", "repetition rate multiplication factor")
    omega: float = param(0.0, "deg", "goniometer rotation about +Y")
    bender_in: bool = param(False, "", "True = bender in, False = bender out")
    b_rot: float = param(BENDER_ROTATION, "deg", "rotation of the cold neutron bender")


DEFAULT_PARAMS = Params()  # frozen, so one shared default is safe


class TRex(Instrument):
    """The T-Rex instrument.

    Parameters
    ----------
    sample:
        The crystal, if any. A copy goes in the instrument, relative to ``sample_pos`` unless
        ``sample.relative`` says otherwise; ``sample`` itself is not changed.
    params:
        Control parameters.
    """

    params_class = Params
    name = "TRex"  # the McStas instrument name

    origin: Progress_bar
    source: ESS_butterfly
    arm_nboa: Arm
    guide_1: Guide_gravity
    guide_2: Guide_gravity
    guide_3: Guide_gravity
    arm_guide_4: Arm
    guide_4a: Guide_gravity
    bender: Pol_bender
    guide_4b: Guide_gravity
    arm_start_of_curve: Arm
    guide_5: Guide_gravity
    guide_6: Guide_gravity
    guide_7: Guide_gravity
    guide_8: Guide_gravity
    guide_9: Guide_gravity
    guide_10: Guide_gravity
    guide_11: Guide_gravity
    guide_12: Guide_gravity
    guide_13a: Guide_gravity
    arm_bwc1: Arm
    guide_13b: Guide_gravity
    guide_14a: Guide_gravity
    guide_14b: Guide_gravity
    guide_15a: Guide_gravity
    arm_bwc2: Arm
    guide_15b: Guide_gravity
    arm_bm1: Arm
    guide_16: Guide_gravity
    guide_17: Guide_gravity
    guide_18: Guide_gravity
    guide_19: Guide_gravity
    guide_20: Guide_gravity
    guide_21: Guide_gravity
    guide_22: Guide_gravity
    guide_23: Guide_gravity
    guide_24: Guide_gravity
    guide_25: Guide_gravity
    guide_26: Guide_gravity
    guide_27: Guide_gravity
    guide_28: Guide_gravity
    guide_29: Guide_gravity
    arm_end_of_curve: Arm
    guide_30: Guide_gravity
    guide_31: Guide_gravity
    guide_32: Guide_gravity
    arm_psc: Arm
    guide_33: Guide_gravity
    guide_34: Guide_gravity
    guide_35: Guide_gravity
    guide_36: Guide_gravity
    guide_37: Guide_gravity
    guide_38: Guide_gravity
    guide_39: Guide_gravity
    guide_40: Guide_gravity
    guide_41: Guide_gravity
    guide_42: Guide_gravity
    guide_43: Guide_gravity
    guide_44: Guide_gravity
    guide_45: Guide_gravity
    guide_47: Guide_gravity
    guide_48b: Guide_gravity
    guide_49: Guide_gravity
    arm_mc: Arm
    guide_50a: Guide_gravity
    slit_mc: Slit
    sample_pos: Arm
    sample: Single_crystal | None

    def __init__(
        self, sample: Single_crystal | None = None, params: Params = DEFAULT_PARAMS
    ) -> None:
        super().__init__(params)
        self.origin = Progress_bar(at=[0, 0, 0])
        self.source = source()
        for guides in (guides_extraction, guides_curve, guides_straight):
            for name, component in guides().items():
                setattr(self, name, component)
        self.slit_mc = slit_mc()
        self.sample_pos = sample_pos()

        self.sample = None
        if sample is not None:
            self.sample = replace(sample, relative=sample.relative or "sample_pos")

    @property
    def instr_params(self) -> list[InstrParam]:
        """McStas instrument parameters: the choice of guide 4 and the rotation of the bender."""
        return [
            InstrParam("bender_in", int(self.params.bender_in), "", "bender in = 1, out = 0"),
            InstrParam("b_rot", self.params.b_rot, "deg", "rotation of the cold neutron bender"),
        ]
