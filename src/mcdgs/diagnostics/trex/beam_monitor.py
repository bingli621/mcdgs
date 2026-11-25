"""The beam at the beam monitor of T-Rex, after guide 4a, from McStas and from tof.

:class:`TRexBeamMonitor` is :class:`.TRex` up to its beam monitor, in the middle of the bunker
gap after guide 4, with the bender out: the beam goes through the thermal guide 4a. The source
and the beam monitor keep their default range, the widest band, ``L_MIN`` to ``L_MAX``, where
those of :class:`.TRex` follow ``lambda0``. It also has the source monitor of
:class:`.TRexSource`, at the entrance of guide 1.

tof has no guides: every neutron of its pulse reaches the beam monitor. Its counts are always
scaled at the source monitor, to the total of McStas there, where both see the same beam. At
the beam monitor, a transmission on top of that scale stands for what the guides take out of
the beam. It is tuned so that tof matches McStas: by default, so that both have the same total.
The transmission is then that of the guides in McStas, and the shape in time and wavelength
that of each model, so what the guides do to the shape shows as a difference between the two.

:func:`.count_rate` reads from the map of McStas the rate a real monitor at that place would
count, on average and at the peak of the pulse, for an efficiency ε(λ) or for the flux itself.

:func:`main` runs both models, reads the intensity of each on the bins of the monitors, scales
tof at the source monitor, tunes the transmission at the beam monitor, and
:func:`.plot_comparison` draws them side by side. It also prints the count rates of the flux and
of a monitor with ε = C λ. Run the module to do it::

    pixi run python -m mcdgs.diagnostics.trex.beam_monitor
"""

from dataclasses import dataclass
from pathlib import Path

from mcdgs.core import Instrument
from mcdgs.core.components import Arm, ESS_butterfly, Guide_gravity, Monitor_nD, Progress_bar
from mcdgs.diagnostics import cli, compare
from mcdgs.diagnostics.beam_monitor import print_rates
from mcdgs.instruments.trex.guides.bispectral import arm_guide_4, guide_4a
from mcdgs.instruments.trex.guides.extraction import guides_extraction
from mcdgs.instruments.trex.monitors import beam_monitor, source_monitor
from mcdgs.instruments.trex.sources import source
from mcdgs.mcstas import Simulation, read_mcstas_maps
from mcdgs.reduction.binning import monitor_bins
from mcdgs.tof import build_tof, read_tof_maps


# --8<-- [start:instrument]
@dataclass(frozen=True)
class TRexBeamMonitorParams:
    """Control parameters of :class:`TRexBeamMonitor`: there are none."""


class TRexBeamMonitor(Instrument):
    """T-Rex from the source to the beam monitor after guide 4, through guide 4a.

    Parameters
    ----------
    params:
        Control parameters. If ``None``, the defaults of :class:`TRexBeamMonitorParams`.
    """

    params_class = TRexBeamMonitorParams
    name = "trex_beam_monitor"  # the McStas instrument name

    origin: Progress_bar
    source: ESS_butterfly
    source_monitor: Monitor_nD
    arm_nboa: Arm
    guide_1: Guide_gravity
    guide_2: Guide_gravity
    guide_3: Guide_gravity
    arm_guide_4: Arm
    guide_4a: Guide_gravity
    beam_monitor: Monitor_nD

    def __init__(self, params: TRexBeamMonitorParams | None = None) -> None:
        super().__init__(params)
        self.origin = Progress_bar(at=[0, 0, 0])
        self.source = source()  # the widest band, as the beam monitor's default limits
        self.source_monitor = source_monitor()
        self.add_components(guides_extraction())
        # the bender out: guide 4a, the thermal guide, always on the beam
        self.arm_guide_4 = arm_guide_4()
        self.guide_4a = guide_4a()
        self.beam_monitor = beam_monitor()


def define_instrument() -> TRexBeamMonitor:
    """Return the instrument of this diagnostic: :class:`TRexBeamMonitor`."""
    params = TRexBeamMonitorParams()  # the default control parameters: there are none
    return TRexBeamMonitor(params)


# --8<-- [end:instrument]


def main(
    out: Path | None = None,
    ncount: float = 1e7,
    neutrons: int = 1_000_000,
    seed: int = 1,
    transmission: float | None = None,
) -> Path:
    """Run both models of :class:`TRexBeamMonitor`, print their means and save the figure.

    The figure is ``trex_beam_monitor.png``.

    Parameters
    ----------
    out:
        Folder of the figure. If ``None``, ``trex_beam_monitor`` in the results folder of
        :meth:`.Workspace.find`.
    ncount:
        Number of McStas neutrons. More than at the source: the guides let through a fraction
        of the beam.
    neutrons:
        Number of tof neutrons.
    seed:
        Seed of the random numbers of both models.
    transmission:
        The factor on tof at the beam monitor, on top of the scale of the source monitor: the
        fraction of the beam the guides let through, which tof does not model. If ``None``,
        tuned so that tof has the total of McStas at the beam monitor.

    Returns
    -------
    :
        The path of the figure.
    """
    instrument = define_instrument()
    name = "trex_beam_monitor"

    # --8<-- [start:run_mcstas]
    # run McStas
    sim = Simulation(instrument)
    mcstas_run = sim.run(folder_name=name, ncount=ncount, seed=seed)
    # --8<-- [end:run_mcstas]

    # --8<-- [start:run_tof]
    # run tof, then read the monitors of both on the bins of the McStas maps: time, wavelength
    tof_model = build_tof(instrument, neutrons=neutrons, seed=seed)
    tof_run = tof_model.run()
    bins = {  # read from the monitors
        name: monitor_bins(getattr(instrument, name).bin_limits)
        for name in ("source_monitor", "beam_monitor")
    }
    mcstas_maps = read_mcstas_maps(mcstas_run, bins)
    tof_maps = read_tof_maps(tof_run, bins)
    # --8<-- [end:run_tof]

    # --8<-- [start:transmission]
    # tof's counts are scaled to the total of McStas at the source monitor, in n/s; at the beam
    # monitor, the same scale, and the transmission on top of it, tuned if it is None
    mcstas_source, tof_source = mcstas_maps["source_monitor"], tof_maps["source_monitor"]
    scale = compare.rescale_tof_to_mcstas(tof_source, mcstas_source)
    at_source = compare.mcstas_vs_tof(
        "source_monitor", mcstas_source, tof_source, ncount, neutrons, scale
    )
    mcstas_beam, tof_beam = mcstas_maps["beam_monitor"], tof_maps["beam_monitor"]
    if transmission is None:  # tuned: the transmission that gives tof the total of McStas
        transmission = compare.rescale_tof_to_mcstas(tof_beam, mcstas_beam) / scale
    comparison = compare.mcstas_vs_tof(
        "beam_monitor", mcstas_beam, tof_beam, ncount, neutrons, scale, transmission
    )
    # --8<-- [end:transmission]

    # print the mean time and wavelength of each model, and the scale and transmission of tof
    compare.print_means(comparison)
    mcstas_total = comparison.first.values.sum()
    print(
        f"McStas total {mcstas_total:.4g} n/s, "
        f"{mcstas_total / at_source.first.values.sum():.4g} of that at the source, "
        f"{scale:.4g} n/s per tof neutron, transmission {transmission:.4g}"
    )
    monitor = instrument.beam_monitor
    print_rates(
        comparison.first.values, bins["beam_monitor"], area=monitor.xwidth * monitor.yheight * 1e4
    )

    fig = compare.plot_comparison(comparison, "T-Rex flux at the beam monitor after guide 4a")
    return cli.save_figure(fig, cli.output_folder(out, name) / f"{name}.png")


if __name__ == "__main__":
    main(out=cli.parse_out(__doc__))
