"""The CSPEC instrument and its control parameters."""

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Any

from mcdgs.core import InstrParam, Instrument, param
from mcdgs.core.components import Arm, Component, Guide_gravity, NXdisk_chopper, Progress_bar
from mcdgs.core.components.masked_source import chopper_train

from .choppers import chopper, chopper_timings
from .guides.bender import guides_bender
from .guides.frame import arm_iscs
from .guides.straight import guides_straight, sample_pos
from .monitors import sample_monitor
from .sources import masked_source, source, wavelength_band


@dataclass(frozen=True)
class Params:
    """Control parameters of CSPEC.

    Frozen: change them with ``CSPEC.update_params``. Each field's unit and comment are in its
    metadata.

    Parameters
    ----------
    t0_offset:
        Time from the start of the pulse at which the neutrons the choppers are phased to leave
        the moderator, in s.
    lambda0:
        Central wavelength, in Å. ``None`` is the white beam: the source simulates the whole
        of ``L_MIN`` to ``L_MAX``. The choppers then have no wavelength to be phased to, so
        only an instrument that parks them can run.
    nu_mc:
        Frequency of the monochromatic choppers, in Hz. The pulse-shaping choppers turn at half
        of it.
    n:
        The RRM chopper turns at ``nu_mc / n``.
    bw3_delay:
        Delay added to that of the bandwidth chopper ``bw3``, in s.
    omega:
        Goniometer rotation about +y, in deg.
    """

    t0_offset: float = param(
        1.7e-3, "s", "time added to choppers shifting to middle of the long pulse"
    )

    lambda0: float | None = param(3.0, "Å", "central wavelength")
    nu_mc: float = param(224.0, "Hz", "frequency of the monochromatic choppers")
    n: int = param(1, "", "the RRM chopper turns at nu_mc / n")
    bw3_delay: float = param(0.0, "s", "delay added to that of chopper bw3")

    omega: float = param(0.0, "deg", "goniometer rotation about +Y")


DEFAULT_PARAMS = Params()  # frozen, so one shared default is safe


class CSPEC(Instrument):
    """The CSPEC instrument.

    Parameters
    ----------
    sample:
        The sample, if any: any component, such as those of :mod:`mcdgs.samples`, relative to
        ``sample_pos`` or to nothing (``relative`` ``None``). A copy goes in the instrument,
        relative to ``sample_pos``; ``sample`` itself is not changed.
    params:
        Control parameters.
    source_masked:
        True uses :func:`.sources.masked_source`, the source masked by the eight choppers, for
        many more neutrons through them per ray; False uses :func:`.sources.source`.

    Raises
    ------
    ValueError
        If ``sample`` is relative to a component other than ``sample_pos``.
    """

    params_class = Params
    name = "CSPEC"  # the McStas instrument name

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
        super().__init__(params)

        to_bwc1, to_bwc2, to_end_of_bender = guides_bender()
        to_bwc3, to_psc, to_mc, to_guide_exchanger = guides_straight()
        # the source simulates the band of lambda0
        band = wavelength_band(self.params.lambda0)

        # the choppers and the masked source need the guides in place: they are set at the end,
        # in the places kept for them here, the choppers by a disk without slits
        self.origin = Progress_bar(at=[0, 0, 0])
        self.source = source(*band)
        # the frame of the sheet, which every guide and gap hangs off
        self.arm_iscs = arm_iscs()
        # the S-bender
        self.add_components(to_bwc1)
        self.bw1 = NXdisk_chopper()
        self.add_components(to_bwc2)
        self.bw2 = NXdisk_chopper()
        self.add_components(to_end_of_bender)
        # the straight guides
        self.add_components(to_bwc3)
        self.bw3 = NXdisk_chopper()
        self.add_components(to_psc)
        self.ps1 = NXdisk_chopper()
        self.ps2 = NXdisk_chopper()
        self.add_components(to_mc)
        self.rrm = NXdisk_chopper()
        self.mc1 = NXdisk_chopper()
        self.mc2 = NXdisk_chopper()
        self.add_components(to_guide_exchanger)
        # sample position, and the beam reaching it
        self.sample_pos = sample_pos()
        self.sample_monitor = sample_monitor(lambda_limits=band)
        self.sample = None if sample is None else replace(sample, relative="sample_pos")

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
            from :func:`.choppers.chopper`, by name, in beam order.
        """
        return {name: chopper(name, guide) for name, guide in self.upstream_guides.items()}

    def update_params(self, **values: Any) -> None:
        """Change control parameters.

        A new ``lambda0`` also gives the source the band of :func:`.sources.wavelength_band`, and
        creates the sample monitor again with it. The source keeps everything else, so a masked
        source stays masked, also in an instrument a parameters file made.

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
            self.sample_monitor = sample_monitor(lambda_limits=band)

    def running_conditions(self) -> str:
        """Describe the settings of CSPEC in one line, such as the title of a figure.

        The name and ``lambda0``, then the signed frequencies of the pulse-shaping (PS), RRM and
        monochromatic (MC) choppers, in beam order.

        Returns
        -------
        :
            Such as ``"CSPEC · λ0 = 3 Å · PS 7/-7 Hz · RRM 14 Hz · MC 14/-14 Hz"``.
        """
        parts = [self.name, f"λ0 = {self.params.lambda0:g} Å"]
        for group in ("ps", "rrm", "mc"):
            frequencies = "/".join(f"{nu:g}" for nu in self.chopper_frequencies(group).values())
            parts.append(f"{group.upper()} {frequencies} Hz")
        return " · ".join(parts)

    @property
    def instr_params(self) -> list[InstrParam]:
        """McStas instrument parameters.

        The rotation of the goniometer, and the frequency and delay of each chopper, from
        :meth:`_chopper_instr_params`.
        """
        gonio_params = [
            InstrParam("omega", self.params.omega, "deg", "rotation of the goniometer about +Y"),
        ]
        return gonio_params + self._chopper_instr_params()

    def _chopper_instr_params(self) -> list[InstrParam]:
        """Return the McStas instrument parameters of the choppers.

        The frequency and delay of each chopper, computed from ``lambda0``, ``nu_mc``, ``n``,
        ``t0_offset`` and ``bw3_delay`` with :func:`.chopper_timings`.

        Returns
        -------
        :
            ``nu_<name>`` and ``delay_<name>`` of each chopper, in beam order.

        Raises
        ------
        ValueError
            If ``lambda0`` is ``None``: the choppers have no wavelength to be phased to.
        """
        if self.params.lambda0 is None:
            raise ValueError("lambda0 is None: the choppers have no wavelength to be phased to")
        timings = chopper_timings(
            self.params.lambda0,
            self.params.nu_mc,
            self.params.n,
            self.source.frequency,
            self.params.t0_offset,
            self.params.bw3_delay,
        )
        return [
            p
            for name, (nu, delay) in timings.items()
            for p in (
                InstrParam(f"nu_{name}", nu, "Hz", f"rotation frequency of chopper {name}"),
                InstrParam(f"delay_{name}", delay, "s", f"delay of chopper {name}"),
            )
        ]
