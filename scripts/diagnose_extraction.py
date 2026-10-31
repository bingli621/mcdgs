"""Create T-Rex up to guide 4, the thermal and the cold path, print it, and write its .instr file
to mcstas_input/trex_extraction/ of the repository.
"""

from dataclasses import dataclass
from pathlib import Path

from mcstasscript.interface.instr import McStas_instr

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
from mcdgs.mcstas.builder import build_mcstas, default_input_dir

WORK_DIR = Path(__file__).resolve().parents[1]  # the repository, wherever the script runs from


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


def build() -> McStas_instr:
    """Create the extraction section, print its summary and return its McStas instrument."""
    extraction = Extraction()
    print(extraction)

    input_path = default_input_dir(extraction, WORK_DIR)
    return build_mcstas(extraction, input_path=input_path)


if __name__ == "__main__":
    trex_extraction = build()
    trex_extraction.write_full_instrument()
    print(f"\nWritten to {trex_extraction.input_path}")

    trex_extraction.show_components()
    trex_extraction.show_diagram()
    # trex_extraction.show_instrument(backend="matplotlib_2d", projection="zx")  # requires matplotlib
