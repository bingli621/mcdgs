"""The T-Rex instrument and its control parameters."""

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Any

from mcdgs.core import InstrParam, Instrument, param
from mcdgs.core.components import (
    Arm,
    Component,
    Guide_gravity,
    NXdisk_chopper,
    Progress_bar,
)
from mcdgs.core.components.masked_source import Polygon_ESS_butterfly, chopper_train

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
        Central wavelength, in Å. ``None`` is the white beam: the source simulates the whole
        of ``L_MIN`` to ``L_MAX``. The choppers then have no wavelength to be phased to, so
        ``rrm`` must be 0, which parks them.
    rrm:
        Repetition rate multiplication factor. 0 parks every chopper, the bandwidth choppers
        too, with the slit of ``mode`` on the beam.
    mode:
        ``"FF"`` puts the wide slits of the pulse-shaping and monochromatic choppers on the beam,
        for high flux; ``"RR"`` the narrow ones, for high resolution; ``"RF"`` the narrow slits
        of the pulse-shaping choppers and the wide ones of the monochromatic choppers, and
        ``"FR"`` the other way around.
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

    lambda0: float | None = param(2.0, "Å", "central wavelength")
    rrm: int = param(18, "", "repetition rate multiplication factor")
    mode: Mode = param(
        "FF", "", "slits of PS then MC: F = wide, high flux, R = narrow, high resolution"
    )

    omega: float = param(0.0, "deg", "goniometer rotation about +Y")
    bender_in: bool = param(False, "", "True = bender in, False = bender out")
    b_rot: float = param(BENDER_ROTATION, "deg", "rotation of the cold neutron bender")


DEFAULT_PARAMS = Params()  # frozen, so one shared default is safe

# the masked source emits where the choppers open as they turn: parked, they do not
MASKED_PARKED = "rrm 0 parks the choppers, and the masked source needs them turning"


class TRex(Instrument):
    """The T-Rex instrument.

    Parameters
    ----------
    sample:
        The sample, if any: any component, such as those of :mod:`mcdgs.samples`, a
        ``Single_crystal`` or a ``Phonon_simple``, relative to ``sample_pos`` or to nothing
        (``relative`` ``None``). A copy goes in the instrument, relative to ``sample_pos``;
        ``sample`` itself is not changed.
    params:
        Control parameters.
    source_masked:
        True uses :func:`.sources.masked_source`, the source masked by the six choppers, for
        many more neutrons through them per ray; False uses :func:`.sources.source`. Not with
        ``rrm`` 0: the masked source needs the choppers turning.

    Raises
    ------
    ValueError
        If ``sample`` is relative to a component other than ``sample_pos``, or if
        ``source_masked`` is True while ``rrm`` is 0.
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
        if sample is not None and sample.relative not in (None, "sample_pos"):
            raise ValueError(f"the sample is relative to {sample.relative!r}, not to 'sample_pos'")
        if source_masked and params.rrm == 0:
            raise ValueError(MASKED_PARKED)
        super().__init__(params)

        to_bwc1, to_bwc2, to_end_of_curve = guides_curve()
        to_psc, to_mc, to_slit = guides_straight()
        # the source simulates the band of lambda0, and the time and wavelength monitors span it
        band = wavelength_band(self.params.lambda0)

        # the choppers and the masked source need the guides in place: they are set at the end,
        # in the places kept for them here, the choppers by a disk without slits
        self.origin = Progress_bar(at=[0, 0, 0])
        self.source = source(*band)
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
        self.bw1 = NXdisk_chopper()
        self.add_components(to_bwc2)
        self.bw2 = NXdisk_chopper()
        self.add_components(to_end_of_curve)
        # straight guides
        self.add_components(to_psc)
        self.ps1 = NXdisk_chopper()
        self.ps2 = NXdisk_chopper()
        self.add_components(to_mc)
        self.mc1 = NXdisk_chopper()
        self.mc2 = NXdisk_chopper()
        self.add_components(to_slit)
        self.slit_mc = slit_mc()
        self.pre_sample_monitor = pre_sample_monitor(lambda_limits=band)
        # before the sample in order, 3 m after it in space: it counts the beam as if there were
        # no sample, and gives the sample each neutron back unchanged
        self.beamstop_monitor = beamstop_monitor(lambda_limits=band)
        # sample position
        self.sample_pos = sample_pos()
        self.sample_monitor = sample_monitor(lambda_limits=band)
        self.sample = None if sample is None else replace(sample, relative="sample_pos")
        # detectors: the beamstop stops the direct beam through the sample
        self.beamstop = beamstop()
        self.banana_detector = banana_detector()

        # each chopper with its flight path, in its place; the masked source is masked by these
        # same choppers
        choppers = self._create_choppers()
        for name, (nx_chopper, _) in choppers.items():
            setattr(self, name, nx_chopper)
        if source_masked:
            self.source = masked_source(chopper_train(list(choppers.values())), *band)

    def _create_choppers(self) -> dict[str, tuple[NXdisk_chopper, float]]:
        """Create the choppers, each with a beam window as wide as the guide in front.

        Returns
        -------
        :
            Each chopper in the place kept for it in :attr:`components` with its flight path,
            from :func:`.choppers.chopper`, by name, in beam order. A disk with no guide in front
            has an unrestricted window.
        """
        return {name: chopper(name, guide) for name, guide in self.upstream_guides.items()}

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
        ValueError
            If ``rrm`` becomes 0 with the masked source; the parameters are then left as they
            were.
        """
        if values.get("rrm") == 0 and isinstance(self.source, Polygon_ESS_butterfly):
            raise ValueError(MASKED_PARKED)
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

        The name, the mode, ``lambda0`` or the white beam, then ``rrm`` and the signed
        frequencies of the pulse-shaping (PS) and monochromatic (MC) choppers, in beam order,
        or that the choppers are parked, with ``rrm`` 0.

        Returns
        -------
        :
            Such as ``"TRex · FF · λ0 = 2 Å · RRM = 18 · PS 189/-189 Hz · MC 252/252 Hz"``, or
            ``"TRex · RR · white beam · choppers parked"``.
        """
        lambda0 = self.params.lambda0
        parts = [
            self.name,
            self.params.mode,
            "white beam" if lambda0 is None else f"λ0 = {lambda0:g} Å",
        ]
        if self.params.rrm == 0:
            parts.append("choppers parked")
            return " · ".join(parts)
        parts.append(f"RRM = {self.params.rrm}")
        for group in ("ps", "mc"):
            frequencies = "/".join(f"{nu:g}" for nu in self.chopper_frequencies(group).values())
            parts.append(f"{group.upper()} {frequencies} Hz")
        return " · ".join(parts)

    @property
    def instr_params(self) -> list[InstrParam]:
        """McStas instrument parameters.

        The choice of guide 4, the rotation of the bender and of the goniometer, and the
        frequency and delay of each chopper, from :meth:`_chopper_instr_params`.
        """
        bender_params = [
            InstrParam("bender_in", int(self.params.bender_in), "", "bender in = 1, out = 0"),
            InstrParam("b_rot", self.params.b_rot, "deg", "rotation of the cold neutron bender"),
        ]
        gonio_params = [
            InstrParam("omega", self.params.omega, "deg", "rotation of the goniometer about +Y"),
        ]

        return bender_params + gonio_params + self._chopper_instr_params()

    def _chopper_instr_params(self) -> list[InstrParam]:
        """Return the McStas instrument parameters of the choppers.

        The frequency, delay and park angle of each chopper, computed from ``lambda0``,
        ``rrm``, ``mode`` and ``t0_offset`` with :func:`.chopper_timings`.

        Returns
        -------
        :
            ``nu_<name>``, ``delay_<name>`` and ``park_<name>`` of each chopper, in beam order.

        Raises
        ------
        ValueError
            If ``lambda0`` is ``None`` while ``rrm`` is not 0: turning choppers have no
            wavelength to be phased to; or if ``rrm`` is 0 with the masked source, as an
            instrument read from a file can have.
        """
        if self.params.rrm == 0 and isinstance(self.source, Polygon_ESS_butterfly):
            raise ValueError(MASKED_PARKED)
        timings = chopper_timings(
            self.params.lambda0,
            self.params.rrm,
            self.source.frequency,
            self.params.mode,
            self.params.t0_offset,
        )
        return [
            p
            for name, timing in timings.items()
            for p in (
                InstrParam(f"nu_{name}", timing.nu, "Hz", f"rotation frequency of chopper {name}"),
                InstrParam(f"delay_{name}", timing.delay, "s", f"delay of chopper {name}"),
                InstrParam(
                    f"park_{name}", timing.park_angle, "deg", f"park angle of chopper {name}"
                ),
            )
        ]
