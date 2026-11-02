"""Create T-Rex up to guide 4, the thermal and the cold path, print it, write its .instr file to
mcstas_input/trex_extraction/ of the repository, and run the diagnostics of McStasScript on the
cold path into mcstas_output/.
"""

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import mcstasscript as ms
from mcstasscript.interface.instr import McStas_instr

from mcdgs.core.base import Instrument
from mcdgs.core.components import Arm, ESS_butterfly, Guide_gravity, Pol_bender, Progress_bar
from mcdgs.core.params import InstrParam, param
from mcdgs.instruments.trex.guides.extraction import BENDER_ROTATION, guides_extraction
from mcdgs.instruments.trex.sources import source
from mcdgs.mcstas.builder import build_mcstas, default_input_dir
from mcdgs.mcstas.run_folders import default_output_dir, run_folder

WORK_DIR = Path(__file__).resolve().parents[2]  # the repository, wherever the script runs from


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
    arm_nboa: Arm
    guide_1: Guide_gravity
    guide_2: Guide_gravity
    guide_3: Guide_gravity
    arm_guide_4: Arm
    guide_4a: Guide_gravity
    bender: Pol_bender
    guide_4b: Guide_gravity

    def __init__(self, params: Params = DEFAULT_PARAMS) -> None:
        super().__init__(params)
        self.origin = Progress_bar(at=[0, 0, 0])
        self.source = source()
        for name, component in guides_extraction().items():
            setattr(self, name, component)

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


def diagnose(instrument: McStas_instr) -> None:
    diag = ms.Diagnostics(instrument)

    output_path = run_folder(default_output_dir(WORK_DIR), None, "diagnose")
    diag.settings(ncount=1e6, suppress_output=True, mpi="auto", output_path=str(output_path))
    diag.show_settings()

    diag.clear_points()
    for component in ("source", "guide_3", "guide_4b"):
        # rays="all": a limited list stops early, which leaves the intensities too low
        diag.add_point(after=component, rays="all")  # pyright: ignore[reportArgumentType]

    diag.run()

    diag.clear_views()
    diag.add_view("x", "y", bins=100)
    # diag.add_view("vx", "vy", bins=200)
    # diag.add_view("dx", "dy", bins=100)
    # diag.add_view("l", bins=100)
    # diag.add_view("t", bins=100)
    diag.add_view("x", "l", bins=100)
    # diag.add_view("t", "l", bins=100)
    diag.add_view("l")
    diag.plot()
    plt.show()  # McStasScript only draws the figure


if __name__ == "__main__":
    trex_extraction = build()
    trex_extraction.write_full_instrument()
    print(f"\nWritten to {trex_extraction.input_path}")

    trex_extraction.show_components()
    trex_extraction.set_parameters(bender_in=1)
    trex_extraction.show_parameters()
    # trex_extraction.show_diagram(variable="l", limits=[0, 6])
    # trex_extraction.show_instrument(backend="matplotlib_2d", projection="zx")

    diagnose(trex_extraction)
