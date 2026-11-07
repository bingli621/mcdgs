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

:func:`main` runs both models, reads the intensity of each on the bins of the monitors, scales
tof at the source monitor, tunes the transmission at the beam monitor, and
:func:`.plot_comparison` draws them side by side. Run the module to do it::

    pixi run python -m mcdgs.diagnostics.trex_04_beam_monitor
"""

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import mcstasscript as ms
import numpy as np
import scipp as sc
from numpy.typing import NDArray

from mcdgs.core.base import Instrument
from mcdgs.core.components import Arm, ESS_butterfly, Guide_gravity, Monitor_nD, Progress_bar
from mcdgs.core.params import InstrParam
from mcdgs.diagnostics.cli import output_folder, parse_out
from mcdgs.diagnostics.compare import plot_comparison, read_comparison
from mcdgs.instruments.trex.guides.bispectral import arm_guide_4, guide_4a
from mcdgs.instruments.trex.guides.extraction import guides_extraction
from mcdgs.instruments.trex.monitors import beam_monitor, source_monitor
from mcdgs.instruments.trex.sources import source
from mcdgs.mcstas.simulation import Simulation
from mcdgs.tof.builder import build_tof
from mcdgs.tof.reader import monitor_bins


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
        for name, component in guides_extraction().items():
            setattr(self, name, component)
        # the bender out: guide 4a, the thermal guide, always on the beam
        self.arm_guide_4 = arm_guide_4()
        self.guide_4a = guide_4a()
        self.beam_monitor = beam_monitor()

    @property
    def instr_params(self) -> list[InstrParam]:
        """McStas instrument parameters: there are none."""
        return []
        # --8<-- [end:trex_beam_monitor]


def _mean(intensity: NDArray, centers: NDArray, axis: int) -> float:
    """Return the mean of a bin center, weighted by the intensity summed over the other axis."""
    return float(np.average(centers, weights=intensity.sum(axis=1 - axis)))


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
    # tof's counts are scaled to the total of McStas at the source monitor, in n/s
    at_source = read_comparison(
        "source_monitor", mcstas_data, tof_result, source_bins, ncount, neutrons
    )
    # --8<-- [end:run_tof]

    # --8<-- [start:transmission]
    # at the beam monitor, the same scale, and the transmission on top of it
    comparison = read_comparison(
        "beam_monitor", mcstas_data, tof_result, bins, ncount, neutrons, scale=at_source.scale
    )
    if transmission is None:  # tuned: the transmission that gives tof the total of McStas
        transmission = float(comparison.mcstas_intensity.sum() / comparison.tof_intensity.sum())
    comparison = read_comparison(
        "beam_monitor",
        mcstas_data,
        tof_result,
        bins,
        ncount,
        neutrons,
        scale=at_source.scale,
        transmission=transmission,
    )
    # --8<-- [end:transmission]

    # print the mean time and wavelength of each model, and the scale and transmission of tof
    time_centers = sc.midpoints(bins["toa"]).to(unit="ms").values
    wavelength_centers = sc.midpoints(bins["wavelength"]).values  # Å
    for label, intensity in (
        ("McStas", comparison.mcstas_intensity),
        ("tof", comparison.tof_intensity),
    ):
        print(
            f"{label:>6}: mean time {_mean(intensity, time_centers, 0):.3f} ms, "
            f"mean wavelength {_mean(intensity, wavelength_centers, 1):.3f} Å"
        )
    mcstas_total = comparison.mcstas_intensity.sum()
    print(
        f"McStas total {mcstas_total:.4g} n/s, "
        f"{mcstas_total / at_source.mcstas_intensity.sum():.4g} of that at the source, "
        f"{comparison.scale:.4g} n/s per tof neutron, transmission {transmission:.4g}"
    )

    figure = output_folder(out, folder) / f"{name}.png"
    fig = plot_comparison(comparison, "T-Rex flux at the beam monitor after guide 4a")
    fig.savefig(figure, dpi=150)
    plt.close(fig)
    print(f"Saved {figure}")
    return figure


if __name__ == "__main__":
    main(out=parse_out(__doc__))
