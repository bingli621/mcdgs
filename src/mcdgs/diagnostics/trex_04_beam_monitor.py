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

:func:`count_rate` reads from the map of McStas the rate a real monitor at that place would
count, on average and at the peak of the pulse, for an efficiency ε(λ) or for the flux itself.

:func:`main` runs both models, reads the intensity of each on the bins of the monitors, scales
tof at the source monitor, tunes the transmission at the beam monitor, and
:func:`.plot_comparison` draws them side by side. It also prints the count rates of the flux and
of a monitor with ε = C λ. Run the module to do it::

    pixi run python -m mcdgs.diagnostics.trex_04_beam_monitor
"""

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path

import mcstasscript as ms
import numpy as np
import scipp as sc
from numpy.typing import NDArray

from mcdgs.core.base import Instrument
from mcdgs.core.components import Arm, ESS_butterfly, Guide_gravity, Monitor_nD, Progress_bar
from mcdgs.core.components.source import ESS_SOURCE_FREQUENCY
from mcdgs.core.params import InstrParam
from mcdgs.diagnostics.cli import output_folder, parse_out, save_figure
from mcdgs.diagnostics.compare import plot_comparison, print_means, read_scaled_comparisons
from mcdgs.instruments.trex.guides.bispectral import arm_guide_4, guide_4a
from mcdgs.instruments.trex.guides.extraction import guides_extraction
from mcdgs.instruments.trex.monitors import beam_monitor, source_monitor
from mcdgs.instruments.trex.sources import source
from mcdgs.mcstas.simulation import Simulation
from mcdgs.tof.builder import build_tof
from mcdgs.tof.reader import monitor_bins

# 1/Å: C of the example efficiency ε = C λ, 1e-4 at 2 Å, within the 1e-5 to 1e-3 of beam
# monitors
EFFICIENCY_PER_ANGSTROM = 5e-5
DEAD_TIME = 100e-9  # s: the dead time of the example monitor


# --8<-- [start:trex_beam_monitor]
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

    @property
    def instr_params(self) -> list[InstrParam]:
        """McStas instrument parameters: there are none."""
        return []
        # --8<-- [end:trex_beam_monitor]


@dataclass(frozen=True)
class CountRate:
    """What a monitor counts, from a time and wavelength map.

    Parameters
    ----------
    average:
        The rate averaged over the pulses, in counts/s.
    peak:
        The highest rate while a pulse passes, in counts/s.
    peak_time:
        When in the pulse it is reached, in ms: the center of its time bin.
    """

    average: float
    peak: float
    peak_time: float


# --8<-- [start:count_rate]
def count_rate(
    intensity: NDArray,
    bins: Mapping[str, sc.Variable],
    efficiency: Callable[[NDArray], NDArray] | None = None,
    frequency: float = ESS_SOURCE_FREQUENCY,
) -> CountRate:
    """Return the average and the peak rate a monitor counts, from a map of the beam on it.

    Each bin of the map is in n/s, averaged over the pulses. Weighted by the efficiency at its
    wavelength, it is in counts/s; summed over wavelength, a time profile. Divided by the
    frequency and by the width of its time bin, a bin of that profile is the rate while the
    pulse passes it.

    Parameters
    ----------
    intensity:
        The map, in n/s per bin, with the shape (time bins, wavelength bins).
    bins:
        Its bin edges, ``"toa"`` and ``"wavelength"``, in that order.
    efficiency:
        The efficiency ε of the monitor at a wavelength in Å, applied to the centers of the
        wavelength bins. If ``None``, 1: the rate of neutrons through the opening.
    frequency:
        The pulse frequency of the source, in Hz.

    Returns
    -------
    :
        The average and the peak rate, and the time of the peak.
    """
    time_edges = bins["toa"].to(unit="s").values
    wavelength_centers = sc.midpoints(bins["wavelength"]).to(unit="Å").values
    counted = intensity if efficiency is None else intensity * efficiency(wavelength_centers)
    profile = counted.sum(axis=1)  # counts/s per time bin, averaged over the pulses
    instantaneous = profile / frequency / np.diff(time_edges)  # counts/s while the pulse passes
    peak_bin = int(np.argmax(instantaneous))
    peak_time = float(sc.midpoints(bins["toa"]).to(unit="ms").values[peak_bin])
    return CountRate(
        average=float(profile.sum()), peak=float(instantaneous[peak_bin]), peak_time=peak_time
    )
    # --8<-- [end:count_rate]


def _print_rates(intensity: NDArray, bins: Mapping[str, sc.Variable], area: float) -> None:
    """Print the rates of the flux, and of the example monitor with ε = C λ."""
    flux = count_rate(intensity, bins)
    print(
        f"flux: {flux.average:.3g} n/s, {flux.average / ESS_SOURCE_FREQUENCY:.3g} n per pulse, "
        f"{flux.average / area:.3g} n/s/cm², peak {flux.peak:.3g} n/s at {flux.peak_time:.2f} ms"
    )
    counted = count_rate(intensity, bins, lambda wavelength: EFFICIENCY_PER_ANGSTROM * wavelength)
    print(
        f"ε = {EFFICIENCY_PER_ANGSTROM:g} λ/Å: {counted.average:.3g} counts/s, "
        f"peak {counted.peak:.3g} counts/s at {counted.peak_time:.2f} ms, "
        f"peak × dead time {counted.peak * DEAD_TIME:.3g}"
    )


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
        Folder of the figure. If ``None``, the folder of the McStas run, next to its
        ``mccode.h5``.
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
    instrument = TRexBeamMonitor()
    name = "trex_beam_monitor"

    # --8<-- [start:run_mcstas]
    # run McStas and load its monitors, in n/s per bin
    folder = Simulation(instrument).run(name, ncount=ncount, seed=seed).parent
    mcstas_data = ms.load_data(str(folder))
    # --8<-- [end:run_mcstas]

    # --8<-- [start:run_tof]
    # run tof, then read the monitors of both on the bins of the McStas maps: time, wavelength
    tof_result = build_tof(instrument, neutrons=neutrons, seed=seed).run()
    source_bins = monitor_bins(instrument.source_monitor)  # read from the monitors
    bins = monitor_bins(instrument.beam_monitor)
    # --8<-- [end:run_tof]

    # --8<-- [start:transmission]
    # tof's counts are scaled to the total of McStas at the source monitor, in n/s; at the beam
    # monitor, the same scale, and the transmission on top of it, tuned if it is None
    comparisons = read_scaled_comparisons(
        mcstas_data,
        tof_result,
        {"source_monitor": source_bins, "beam_monitor": bins},
        ncount,
        neutrons,
        transmission=transmission,
    )
    at_source, comparison = comparisons["source_monitor"], comparisons["beam_monitor"]
    # --8<-- [end:transmission]

    # print the mean time and wavelength of each model, and the scale and transmission of tof
    print_means(comparison)
    mcstas_total = comparison.mcstas_intensity.sum()
    print(
        f"McStas total {mcstas_total:.4g} n/s, "
        f"{mcstas_total / at_source.mcstas_intensity.sum():.4g} of that at the source, "
        f"{comparison.scale:.4g} n/s per tof neutron, transmission {comparison.transmission:.4g}"
    )
    monitor = instrument.beam_monitor
    _print_rates(comparison.mcstas_intensity, bins, area=monitor.xwidth * monitor.yheight * 1e4)

    fig = plot_comparison(comparison, "T-Rex flux at the beam monitor after guide 4a")
    return save_figure(fig, output_folder(out, folder) / f"{name}.png")


if __name__ == "__main__":
    main(out=parse_out(__doc__))
