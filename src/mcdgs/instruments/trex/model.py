"""The T-Rex instrument and its control parameters."""

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Any

from mcdgs.core.base import Instrument
from mcdgs.core.components import (
    Arm,
    Component,
    Guide_gravity,
    Progress_bar,
)
from mcdgs.core.params import InstrParam, param

from .choppers import Mode, chopper, chopper_timings
from .detectors import banana_detector
from .guides.bispectral import (
    BENDER_ROTATION,
    arm_guide_4,
    bender,
    guide_4a,
    guide_4b,
)
from .guides.curve import guides_curve
from .guides.extraction import guides_extraction
from .guides.straight import guides_straight, sample_pos, slit_mc
from .monitors import (
    beam_monitor,
    beamstop,
    beamstop_monitor,
    pre_sample_monitor,
    sample_monitor,
    source_monitor,
)
from .sources import masked_source, source, wavelength_band


@dataclass(frozen=True)
class Params:
    """Control parameters of T-Rex.

    Frozen: change them with ``TRex.update_params``. Each field's unit and comment are in its
    metadata.

    Parameters
    ----------
    t0_offset:
        Time from the start of the pulse at which the neutrons the choppers are phased to leave
        the moderator, in s.
    lambda0:
        Central wavelength, in Å.
    rrm:
        Repetition rate multiplication factor.
    mode:
        ``"FF"`` puts the wide slits of the pulse-shaping and monochromatic choppers on the beam,
        for high flux; ``"RR"`` the narrow ones, for high resolution.
    omega:
        Goniometer rotation about +y, in deg.
    bender_in:
        True uses the bender and the cold guide 4b, False uses the thermal guide 4a.
    b_rot:
        Rotation of the bender about +y, in deg.
    """

    t0_offset: float = param(
        0.0017, "s", "time added to choppers shifting to middle of the long pulse"
    )

    lambda0: float = param(2.0, "Å", "central wavelength")
    rrm: int = param(18, "", "repetition rate multiplication factor")
    mode: Mode = param("FF", "", "FF = high flux, RR = high resolution")

    omega: float = param(0.0, "deg", "goniometer rotation about +Y")
    bender_in: bool = param(False, "", "True = bender in, False = bender out")
    b_rot: float = param(BENDER_ROTATION, "deg", "rotation of the cold neutron bender")


DEFAULT_PARAMS = Params()  # frozen, so one shared default is safe


class TRex(Instrument):
    """The T-Rex instrument.

    Parameters
    ----------
    sample:
        The sample, if any: any component, such as those of :mod:`mcdgs.samples`, a
        ``Single_crystal`` or a ``Phonon_simple``. A copy goes in the instrument, relative to
        ``sample_pos`` unless ``sample.relative`` says otherwise; ``sample`` itself is not
        changed.
    params:
        Control parameters.
    source_masked:
        True uses :func:`.sources.masked_source`, the source masked by the six choppers, for
        many more neutrons through them per ray; False uses :func:`.sources.source`.
    """

    params_class = Params
    name = "TRex"  # the McStas instrument name

    if TYPE_CHECKING:  # the arms and the guides are set by name: see add_components

        def __getattr__(self, name: str) -> Arm | Guide_gravity:
            """Type the attributes that ``__init__`` sets by name."""
            ...

    def __init__(
        self,
        sample: Component | None = None,
        params: Params = DEFAULT_PARAMS,
        source_masked: bool = False,
    ) -> None:
        super().__init__(params)

        to_bwc1, to_bwc2, to_end_of_curve = guides_curve()
        to_psc, to_mc, to_slit = guides_straight()

        # the source simulates the band of lambda0, and the time and wavelength monitors span it
        band = wavelength_band(self.params.lambda0)

        self.origin = Progress_bar(at=[0, 0, 0])
        self.source = (masked_source if source_masked else source)(*band)
        self.source_monitor = source_monitor()  # the entrance of guide 1, the widest band
        self.add_components(guides_extraction())
        # bispectral switch: guide 4a (thermal) or bender and guide 4b (cold)
        self.arm_guide_4 = arm_guide_4()
        self.guide_4a = guide_4a(when="bender_in==0")
        self.bender = bender(when="bender_in==1")
        self.guide_4b = guide_4b(when="bender_in==1")
        self.beam_monitor = beam_monitor(lambda_limits=band)
        # curved guides
        self.add_components(to_bwc1)
        self.bw1 = chopper("bw1")
        self.add_components(to_bwc2)
        self.bw2 = chopper("bw2")
        self.add_components(to_end_of_curve)
        # straight guides
        self.add_components(to_psc)
        self.ps1 = chopper("ps1")
        self.ps2 = chopper("ps2")
        self.add_components(to_mc)
        self.mc1 = chopper("mc1")
        self.mc2 = chopper("mc2")
        self.add_components(to_slit)
        self.slit_mc = slit_mc()
        self.pre_sample_monitor = pre_sample_monitor(lambda_limits=band)
        # before the sample in order, 3 m after it in space: it counts the beam as if there were
        # no sample, and gives the sample each neutron back unchanged
        self.beamstop_monitor = beamstop_monitor(lambda_limits=band)
        # sample position
        self.sample_pos = sample_pos()
        self.sample_monitor = sample_monitor(lambda_limits=band)
        self.sample = None
        if sample is not None:
            self.sample = replace(sample, relative=sample.relative or "sample_pos")
        # detectors: the beamstop stops the direct beam through the sample
        self.beamstop = beamstop()
        self.banana_detector = banana_detector()

    def update_params(self, **values: Any) -> None:
        """Change control parameters.

        A new ``lambda0`` also gives the source the band of :func:`.sources.wavelength_band`,
        and creates the time and wavelength monitors again (the beam, pre-sample, beamstop and
        sample monitors) with it. The source keeps everything else, so a masked source stays
        masked, also in an instrument a parameters file made.

        Parameters
        ----------
        **values:
            New values by parameter name, e.g. ``update_params(lambda0=4.0)``.

        Raises
        ------
        TypeError
            If a name is not a field of :class:`Params`.
        """
        super().update_params(**values)
        if "lambda0" in values:
            band = wavelength_band(self.params.lambda0)
            self.source = replace(self.source, Lmin=band[0], Lmax=band[1])
            self.beam_monitor = beam_monitor(lambda_limits=band)
            self.pre_sample_monitor = pre_sample_monitor(lambda_limits=band)
            self.beamstop_monitor = beamstop_monitor(lambda_limits=band)
            self.sample_monitor = sample_monitor(lambda_limits=band)

    def running_conditions(self) -> str:
        """Describe the settings of T-Rex in one line, such as the title of a figure.

        The name, the mode, ``lambda0`` and ``rrm``, then the signed frequencies of the
        pulse-shaping (PS) and monochromatic (MC) choppers, in beam order.

        Returns
        -------
        :
            Such as ``"TRex · FF · λ0 = 2 Å · RRM = 18 · PS 189/-189 Hz · MC 252/252 Hz"``.
        """
        parts = [
            self.name,
            self.params.mode,
            f"λ0 = {self.params.lambda0:g} Å",
            f"RRM = {self.params.rrm}",
        ]
        for group in ("ps", "mc"):
            frequencies = "/".join(f"{nu:g}" for nu in self.chopper_frequencies(group).values())
            parts.append(f"{group.upper()} {frequencies} Hz")
        return " · ".join(parts)

    @property
    def instr_params(self) -> list[InstrParam]:
        """McStas instrument parameters.

        The choice of guide 4, the rotation of the bender and of the goniometer, and the
        frequency and delay of each chopper, computed from ``lambda0``, ``rrm``, ``mode`` and
        ``t0_offset``.
        """
        timings = chopper_timings(
            self.params.lambda0,
            self.params.rrm,
            self.source.frequency,
            self.params.mode,
            self.params.t0_offset,
        )
        bender_params = [
            InstrParam("bender_in", int(self.params.bender_in), "", "bender in = 1, out = 0"),
            InstrParam("b_rot", self.params.b_rot, "deg", "rotation of the cold neutron bender"),
        ]

        chopper_params = [
            p
            for name, timing in timings.items()
            for p in (
                InstrParam(f"nu_{name}", timing.nu, "Hz", f"rotation frequency of chopper {name}"),
                InstrParam(f"delay_{name}", timing.delay, "s", f"delay of chopper {name}"),
            )
        ]
        gonio_params = [
            InstrParam("omega", self.params.omega, "deg", "rotation of the goniometer about +Y"),
        ]

        return bender_params + gonio_params + chopper_params
