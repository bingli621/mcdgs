"""The T-Rex instrument and its control parameters."""

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING

from mcdgs.core.base import Instrument
from mcdgs.core.components import (
    Arm,
    Component,
    Guide_gravity,
    Progress_bar,
    Single_crystal,
)
from mcdgs.core.params import InstrParam, param

from .choppers import choppers
from .guides.bispectral import BENDER_ROTATION, arm_guide_4, bender, guide_4a, guide_4b
from .guides.curve import guides_curve
from .guides.extraction import guides_extraction
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

    if TYPE_CHECKING:  # the arms and the guides are set by name: see add_components

        def __getattr__(self, name: str) -> Arm | Guide_gravity:
            """Type the attributes that ``__init__`` sets by name."""
            ...

    def __init__(
        self, sample: Single_crystal | None = None, params: Params = DEFAULT_PARAMS
    ) -> None:
        super().__init__(params)

        disks = choppers()
        to_bwc1, to_bwc2, to_end_of_curve = guides_curve()
        to_psc, to_mc, to_slit = guides_straight()

        self.origin = Progress_bar(at=[0, 0, 0])
        self.source = source()
        self.add_components(guides_extraction())
        # bispectral switch: guide 4a (thermal) or bender and guide 4b (cold)
        self.arm_guide_4 = arm_guide_4()
        self.guide_4a = guide_4a(when="bender_in==0")
        self.bender = bender(when="bender_in==1")
        self.guide_4b = guide_4b(when="bender_in==1")
        # curved guides
        self.add_components(to_bwc1)
        self.bw1 = disks["bw1"]
        self.add_components(to_bwc2)
        self.bw2 = disks["bw2"]
        self.add_components(to_end_of_curve)
        # straight guides
        self.add_components(to_psc)
        self.ps1 = disks["ps1"]
        self.ps2 = disks["ps2"]
        self.add_components(to_mc)
        self.mc1 = disks["mc1"]
        self.mc2 = disks["mc2"]
        self.add_components(to_slit)
        self.slit_mc = slit_mc()
        #  sample position
        self.sample_pos = sample_pos()
        self.sample = None
        if sample is not None:
            self.sample = replace(sample, relative=sample.relative or "sample_pos")

    def add_components(self, components: dict[str, Component]) -> None:
        """Set components as attributes under their names, keeping their order.

        Parameters
        ----------
        components:
            Components by name, in beam order.
        """
        for name, component in components.items():
            setattr(self, name, component)

    @property
    def instr_params(self) -> list[InstrParam]:
        """McStas instrument parameters: the choice of guide 4 and the rotation of the bender."""
        return [
            InstrParam("bender_in", int(self.params.bender_in), "", "bender in = 1, out = 0"),
            InstrParam("b_rot", self.params.b_rot, "deg", "rotation of the cold neutron bender"),
            InstrParam("omega", self.params.omega, "deg", "rotation of the goniometer about +Y"),
            # InstrParam("nu_bw1", self.params.nu_bw1, "Hz", "frequency of bandwidth chopper 1"),
            # InstrParam("nu_bw2", self.params.nu_bw2, "Hz", "frequency of bandwidth chopper 2"),
        ]
