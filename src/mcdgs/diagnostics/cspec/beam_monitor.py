"""The beam at the beam monitor of CSPEC, after the BBG, from McStas and from tof.

:class:`CSPECBeamMonitor` is :class:`.CSPEC` up to its beam monitor, in the middle of the gap
after the BBG, row 22, 5.9 m from the moderator: the source, the source monitor at the entrance
of W03-01-011, and the guides from W03-01-011 to the BBG, W03-02-01. No chopper stands before
it. The source and the beam monitor keep the widest band, ``L_MIN`` to ``L_MAX``, where those of
:class:`.CSPEC` follow ``lambda0``.

As at the beam monitor of T-Rex, tof has no guides: every neutron of its pulse reaches the beam
monitor. Its counts are scaled at the source monitor, to the total of McStas there, where both
see the same beam. At the beam monitor, a transmission on top of that scale stands for what the
guides take out of the beam. It is tuned by default so that both have the same intensity at
3 Å, the default ``lambda0`` of :class:`.CSPEC`: summed over time, in the wavelength bins
within 0.05 Å of it (:func:`intensity_near`). Away from 3 Å, a difference between the two is
what the guides do to the beam, relative to what they do at 3 Å.

:func:`.count_rate` reads from the map of McStas the rate a real monitor at that place would
count, on average and at the peak of the pulse.

:func:`main` runs both models, reads the intensity of each on the bins of the monitors, scales
tof at the source monitor, tunes the transmission at the beam monitor, and
:func:`.plot_comparison` draws them side by side. It also prints the count rates of the flux and
of a monitor with ε = C λ. Run the module to do it::

    pixi run python -m mcdgs.diagnostics.cspec.beam_monitor
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipp as sc

from mcdgs.core import Instrument
from mcdgs.core.components import Arm, ESS_butterfly, Guide_gravity, Monitor_nD, Progress_bar
from mcdgs.diagnostics import cli, compare
from mcdgs.diagnostics.beam_monitor import print_rates
from mcdgs.instruments.cspec.guides.bender import guides_bender
from mcdgs.instruments.cspec.guides.frame import arm_iscs
from mcdgs.instruments.cspec.monitors import beam_monitor, source_monitor
from mcdgs.instruments.cspec.sources import source
from mcdgs.mcstas import Simulation, read_mcstas_maps
from mcdgs.reduction.binning import monitor_bins
from mcdgs.tof import build_tof, read_tof_maps

# Å: the wavelength at which tof is matched to McStas, the default lambda0 of CSPEC, and the
# half width of the window around it
MATCH_WAVELENGTH = 3.0
MATCH_HALF_WIDTH = 0.05


# --8<-- [start:instrument]
@dataclass(frozen=True)
class CSPECBeamMonitorParams:
    """Control parameters of :class:`CSPECBeamMonitor`: there are none."""


class CSPECBeamMonitor(Instrument):
    """CSPEC from the source to the beam monitor in the gap after the BBG.

    Parameters
    ----------
    params:
        Control parameters. If ``None``, the defaults of :class:`CSPECBeamMonitorParams`.
    """

    params_class = CSPECBeamMonitorParams
    name = "cspec_beam_monitor"  # the McStas instrument name

    origin: Progress_bar
    source: ESS_butterfly
    arm_iscs: Arm
    source_monitor: Monitor_nD
    guide_01_011: Guide_gravity
    guide_01_012: Guide_gravity
    guide_01_02: Guide_gravity
    guide_01_03: Guide_gravity
    guide_02_01: Guide_gravity
    arm_gap_bbg: Arm
    beam_monitor: Monitor_nD

    def __init__(self, params: CSPECBeamMonitorParams | None = None) -> None:
        super().__init__(params)
        to_gap_bbg, *_ = guides_bender()
        self.origin = Progress_bar(at=[0, 0, 0])
        self.source = source()  # the widest band, as the beam monitor's default limits
        self.arm_iscs = arm_iscs()  # the frame of the sheet, which every guide hangs off
        self.source_monitor = source_monitor()
        # W03-01-011 to the BBG, and the arm at the start of the gap after it
        self.add_components(to_gap_bbg)
        self.beam_monitor = beam_monitor()


def define_instrument() -> CSPECBeamMonitor:
    """Return the instrument of this diagnostic: :class:`CSPECBeamMonitor`."""
    params = CSPECBeamMonitorParams()  # the default control parameters: there are none
    return CSPECBeamMonitor(params)


# --8<-- [end:instrument]


def intensity_near(
    intensity: sc.DataArray, wavelength: float, half_width: float = MATCH_HALF_WIDTH
) -> float:
    """Return the intensity of a map summed over time, in the bins near a wavelength.

    Parameters
    ----------
    intensity:
        A time and wavelength map, with the bin edges ``"wavelength"``.
    wavelength:
        The wavelength, in Å.
    half_width:
        The bins whose centers are within this of ``wavelength``, in Å, are summed.

    Returns
    -------
    :
        The sum, in the unit of the map.

    Raises
    ------
    ValueError
        If no bin center is within ``half_width`` of ``wavelength``.
    """
    centers = sc.midpoints(intensity.coords["wavelength"]).to(unit="Å").values
    (near,) = np.nonzero(np.abs(centers - wavelength) <= half_width)
    if near.size == 0:
        raise ValueError(f"no wavelength bin within {half_width} Å of {wavelength} Å")
    return float(intensity["wavelength", near[0] : near[-1] + 1].sum().value)


def main(
    out: Path | None = None,
    ncount: float = 1e7,
    neutrons: int = 1_000_000,
    seed: int = 1,
    transmission: float | None = None,
    match_wavelength: float = MATCH_WAVELENGTH,
) -> Path:
    """Run both models of :class:`CSPECBeamMonitor`, print their means and save the figure.

    The figure is ``cspec_beam_monitor.png``.

    Parameters
    ----------
    out:
        Folder of the figure. If ``None``, ``cspec_beam_monitor`` in the results folder of
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
        tuned so that tof has the intensity of McStas at ``match_wavelength``.
    match_wavelength:
        The wavelength, in Å, at which ``transmission`` is tuned: see :func:`intensity_near`.

    Returns
    -------
    :
        The path of the figure.
    """
    instrument = define_instrument()
    name = "cspec_beam_monitor"

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
    if transmission is None:  # tuned: tof gets the intensity of McStas at match_wavelength
        mcstas_near = intensity_near(mcstas_beam, match_wavelength)
        tof_near = intensity_near(tof_beam, match_wavelength)
        transmission = mcstas_near / (tof_near * scale)
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
        f"tof total {comparison.second.values.sum():.4g} n/s, "
        f"{scale:.4g} n/s per tof neutron, transmission {transmission:.4g}"
    )
    monitor = instrument.beam_monitor
    print_rates(
        comparison.first.values, bins["beam_monitor"], area=monitor.xwidth * monitor.yheight * 1e4
    )

    fig = compare.plot_comparison(comparison, "CSPEC flux at the beam monitor after the BBG")
    return cli.save_figure(fig, cli.output_folder(out, name) / f"{name}.png")


if __name__ == "__main__":
    main(out=cli.parse_out(__doc__))
