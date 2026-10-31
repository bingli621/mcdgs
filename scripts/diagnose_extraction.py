"""Create T-Rex up to guide 4, the thermal and the cold path, and print it."""

from dataclasses import dataclass

from mcdgs.core.base import Instrument
from mcdgs.core.components import (
    Arm,
    ESS_butterfly,
    Guide_gravity,
    Monitor_nD,
    Pol_bender,
    Progress_bar,
)
from mcdgs.core.params import InstrParam, param
from mcdgs.instruments.trex.guides_extraction import BENDER_ROTATION, guides_extraction
from mcdgs.instruments.trex.trex import source


@dataclass(frozen=True)
class Params:
    """Control parameters of the extraction section.

    Frozen: change them with ``Extraction.update_params``. Each field's unit and comment are in
    its metadata.

    Parameters
    ----------
    bender_in:
        True uses the bender and the cold guide 4b, False uses the thermal guide 4a.
    b_rot:
        Rotation of the bender about +y, in deg.

    """

    bender_in: bool = param(False, "", "True = bender in, False = bender out")
    b_rot: float = param(BENDER_ROTATION, "deg", "rotation of the cold neutron bender")


DEFAULT_PARAMS = Params()  # frozen, so one shared default is safe


class Extraction(Instrument):
    """T-Rex up to guide 4, the thermal and the cold path.

    Parameters
    ----------
    params:
        Control parameters.
    """

    params_class = Params
    name = "trex_extraction"  # the McStas instrument name

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
    diag_after_guide_4b: Monitor_nD

    def __init__(self, params: Params = DEFAULT_PARAMS) -> None:
        super().__init__(params)
        self.origin = Progress_bar(at=[0, 0, 0])
        self.source = source()
        for name, component in guides_extraction().items():
            setattr(self, name, component)
        self.diag_after_guide_4b = Monitor_nD(  # as the Diag monitors of trex_primary.txt
            xwidth=100,
            yheight=100,
            restore_neutron=1,
            options="previous n x y z vx vy vz t, list 50001",
            filename="diag_after_guide_4b.diag",
            at=[0, 0, self.guide_4b.l + 1e-3],  # the shared exit of guides 4a and 4b
            relative="guide_4b",
        )

    @property
    def instr_params(self) -> list[InstrParam]:
        return [
            InstrParam("bender_in", int(self.params.bender_in), "", "bender in = 1, out = 0"),
            InstrParam("b_rot", self.params.b_rot, "deg", "rotation of the cold neutron bender"),
        ]


def main() -> None:
    """Create the extraction section and print its summary."""
    print(Extraction())


if __name__ == "__main__":
    main()
