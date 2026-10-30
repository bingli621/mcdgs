from dataclasses import dataclass, replace

from mcdgs.core.base import Instrument
from mcdgs.core.components.crystal import Single_crystal
from mcdgs.core.components.misc import Arm, Progress_bar
from mcdgs.core.components.optics import Guide_gravity, Pol_bender
from mcdgs.core.components.source import ESS_butterfly
from mcdgs.core.params import InstrParam, param

from .guides_extraction import BENDER_ROTATION, guides_extraction

SAMPLE_X = 0.8301227276790154  # m, current estimate, not in the sheet
SAMPLE_Z = 163.8  # m, current estimate, cell C148
SAMPLE_ROTATION = 0.430122381481144  # deg, current estimate, not in the sheet


@dataclass(frozen=True)
class Params:
    """Control parameters of the diffractometer.

    Frozen: change them with ``SingleCrystalDiffractometer.update_params``. Each field's unit and
    comment are in its metadata.

    Parameters
    ----------
    omega:
        Goniometer rotation about +y, in deg.
    bender:
        0 uses the thermal guide 4a, 1 the bender and the cold guide 4b.
    b_rot:
        Rotation of the bender about +y, in deg.

    """

    lambda0: float = param(2.0, "Å", "central wavelength")
    rrm: int = param(18, "", "repetition rate multiplication factor")
    omega: float = param(0.0, "deg", "goniometer rotation about +Y")
    bender: bool = param(False, "", "bender out = 0, bender in = 1")
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
    name = "T-Rex"  # the McStas instrument name

    origin: Progress_bar
    source: ESS_butterfly
    arm_NBOA: Arm
    guide_1: Guide_gravity
    guide_2: Guide_gravity
    guide_3: Guide_gravity
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
        self.source = ESS_butterfly(  # as trex_primary.txt
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
        for name, component in guides_extraction().items():
            setattr(self, name, component)

        self.sample_pos = Arm(at=[SAMPLE_X, 0, SAMPLE_Z], rotated=[0, SAMPLE_ROTATION, 0])
        self.sample = None
        if sample is not None:
            self.sample = replace(sample, relative=sample.relative or "sample_pos")

    @property
    def instr_params(self) -> list[InstrParam]:
        return [
            InstrParam("bender", int(self.params.bender), "", "bender out = 0, bender in = 1"),
            InstrParam("b_rot", self.params.b_rot, "deg", "rotation of the cold neutron bender"),
        ]
