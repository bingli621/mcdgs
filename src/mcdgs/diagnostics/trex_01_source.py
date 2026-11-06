"""The beam at the entrance of guide 1 of T-Rex, from McStas and from tof.

:class:`TRexSource` is the source and a time and wavelength monitor at the entrance of
guide 1, the first optics of T-Rex, 1.895927 m from the source with the opening of that guide.
Nothing stands between the two, so McStas and tof should see the same pulse there: only their
models of the ESS moderator differ.

Counts and intensity. Each simulated neutron stands for some of the real flux, its weight:
intensity = counts × weight, in n/s. The counts are how many simulated neutrons land in a bin,
and depend only on how many are simulated; the intensity is the physics, and should not.

- McStas sets the weights itself, so that the rays add up to the flux of the moderator. With
  ten times the ``ncount``, a bin has ten times the rays, each of a tenth of the weight: the
  intensity stays, only its noise goes down.
- tof's pulse is a distribution with no flux: each neutron has weight 1, so its intensity is
  its count, and grows with the number of neutrons drawn. So tof takes its weight from McStas:
  ``scale``, the total of McStas over that of tof, is the n/s each tof neutron stands for.
  Scaled, the level is that of McStas, and the shape in time and wavelength that of each
  model.

Error bars. A bin of N neutrons has the Poisson error √N, a relative error of 1/√N. In n/s it
is √N × weight. A hundred times the neutrons gives ten times the error in counts, but a tenth
of the relative error, and a tenth of the error in n/s, as the weight is a hundredth. So the
error bars shrink with more neutrons, in both models, while the intensity stays.

:func:`main` runs both models, reads the intensity of each on the bins of the monitor, scales
tof to McStas, and :func:`.plot_comparison` draws them side by side. Run the module to do all
of it::

    pixi run python -m mcdgs.diagnostics.trex_01_source
"""

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import mcstasscript as ms
import numpy as np
import scipp as sc
from numpy.typing import NDArray

from mcdgs.core.base import Instrument
from mcdgs.core.components import ESS_butterfly, Monitor_nD
from mcdgs.core.params import InstrParam
from mcdgs.diagnostics.cli import output_folder, parse_out
from mcdgs.diagnostics.compare import plot_comparison, read_comparison
from mcdgs.instruments.trex.monitors import source_monitor
from mcdgs.instruments.trex.sources import source
from mcdgs.mcstas.simulation import Simulation
from mcdgs.tof.builder import build_tof
from mcdgs.tof.reader import monitor_bins


# --8<-- [start:trex_source]
@dataclass(frozen=True)
class TRexSourceParams:
    """Control parameters of :class:`TRexSource`: there are none."""


class TRexSource(Instrument):
    """The source of T-Rex and a time and wavelength monitor at the entrance of guide 1.

    The monitor is the entrance of guide 1: its opening, where guide 1 begins, 2.0 mm off the axis
    of the source and turned with the guide.

    Parameters
    ----------
    params:
        Control parameters. If ``None``, the defaults of :class:`TRexSourceParams`.
    """

    params_class = TRexSourceParams
    name = "trex_source"  # the McStas instrument name

    source: ESS_butterfly
    source_monitor: Monitor_nD

    def __init__(self, params: TRexSourceParams | None = None) -> None:
        super().__init__(params)
        self.source = source()
        self.source_monitor = source_monitor()

    @property
    def instr_params(self) -> list[InstrParam]:
        """McStas instrument parameters: there are none."""
        return []
        # --8<-- [end:trex_source]


def _mean(intensity: NDArray, centers: NDArray, axis: int) -> float:
    """Return the mean of a bin center, weighted by the intensity summed over the other axis."""
    return float(np.average(centers, weights=intensity.sum(axis=1 - axis)))


def main(
    out: Path | None = None,
    ncount: float = 1e6,
    neutrons: int = 1_000_000,
    seed: int = 1,
) -> Path:
    """Run both models of :class:`TRexSource`, print their means and save the figure.

    Parameters
    ----------
    out:
        Folder of the figure, ``trex_source.png``. If ``None``, the folder of the McStas run,
        next to its ``mccode.h5``.
    ncount:
        Number of McStas neutrons.
    neutrons:
        Number of tof neutrons.
    seed:
        Seed of the random numbers of both models.

    Returns
    -------
    :
        The path of the figure.
    """
    instrument = TRexSource()

    # --8<-- [start:run_mcstas]
    # run McStas and load its monitors, in n/s per bin
    folder = Simulation(instrument).run("trex_source", ncount=ncount, seed=seed).parent
    mcstas_data = ms.load_data(str(folder))
    # --8<-- [end:run_mcstas]

    # --8<-- [start:run_tof]
    # run tof, then read the monitor of both on the bins of the McStas map: time, wavelength
    tof_result = build_tof(instrument, neutrons=neutrons, seed=seed).run()
    bins = monitor_bins(instrument.source_monitor)  # read from the monitor
    # tof's counts are scaled to the total of McStas, in n/s
    comparison = read_comparison("source_monitor", mcstas_data, tof_result, bins, ncount, neutrons)
    # --8<-- [end:run_tof]

    # print the mean time and wavelength of each model, and the total tof is scaled to
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
    print(
        f"McStas total {comparison.mcstas_intensity.sum():.4g} n/s, "
        f"{comparison.scale:.4g} n/s per tof neutron"
    )

    figure = output_folder(out, folder) / "trex_source.png"
    fig = plot_comparison(comparison, "T-Rex flux at the entrance of guide 1")
    fig.savefig(figure, dpi=150)
    plt.close(fig)
    print(f"Saved {figure}")
    return figure


if __name__ == "__main__":
    main(out=parse_out(__doc__))
