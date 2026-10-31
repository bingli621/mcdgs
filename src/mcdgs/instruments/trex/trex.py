"""The T-Rex instrument, its control parameters and its source."""

from dataclasses import dataclass, replace

from mcdgs.core.base import Instrument
from mcdgs.core.components import (
    Arm,
    ESS_butterfly,
    Guide_gravity,
    Pol_bender,
    Progress_bar,
    Single_crystal,
)
from mcdgs.core.params import InstrParam, param

from .guides_extraction import BENDER_ROTATION, guides_extraction

SAMPLE_X = 0.8301227276790154  # m, current estimate, not in the sheet
SAMPLE_Z = 163.8  # m, current estimate, cell C148
SAMPLE_ROTATION = 0.430122381481144  # deg, current estimate, not in the sheet


@dataclass(frozen=True)
class Params:
    """Control parameters of the diffractometer.

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


def source() -> ESS_butterfly:
    """Create the source component for T-Rex.

    Returns
    -------
    :
        A new source.
    """
    return ESS_butterfly(
        sector="W",
        beamline=7,
        yheight=0.03,
        cold_frac=0.3,
        dist=2,
        focus_xw=0.095,
        focus_yh=0.035,
        Lmin=0.1,
        Lmax=10.0,
        n_pulses=1,
        acc_power=2,
    )


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
    arm_NBOA: Arm
    guide_1: Guide_gravity
    guide_2: Guide_gravity
    guide_3: Guide_gravity
    arm_guide_4: Arm
    guide_4a: Guide_gravity
    bender: Pol_bender
    guide_4b: Guide_gravity
    sample_pos: Arm
    sample: Single_crystal | None

    def __init__(
        self, sample: Single_crystal | None = None, params: Params = DEFAULT_PARAMS
    ) -> None:
        super().__init__(params)
        self.origin = Progress_bar(at=[0, 0, 0])
        self.source = source()
        for name, component in guides_extraction().items():
            setattr(self, name, component)

        self.sample_pos = Arm(at=[SAMPLE_X, 0, SAMPLE_Z], rotated=[0, SAMPLE_ROTATION, 0])
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
