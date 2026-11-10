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
of it, at 1e6 and at 1e5 neutrons in both models, and with the band of ``lambda0`` = 2 Å::

    pixi run python -m mcdgs.diagnostics.trex_01_source
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import mcstasscript as ms

from mcdgs.core.base import Instrument
from mcdgs.core.components import ESS_butterfly, Monitor_nD
from mcdgs.core.params import InstrParam, param
from mcdgs.diagnostics.cli import output_folder, parse_out, save_figure
from mcdgs.diagnostics.compare import plot_comparison, print_means, read_comparison
from mcdgs.instruments.trex.monitors import source_monitor
from mcdgs.instruments.trex.sources import source, wavelength_band
from mcdgs.mcstas.simulation import Simulation
from mcdgs.tof.builder import build_tof
from mcdgs.tof.reader import monitor_bins


# --8<-- [start:trex_source]
@dataclass(frozen=True)
class TRexSourceParams:
    """Control parameters of :class:`TRexSource`.

    Parameters
    ----------
    lambda0:
        Central wavelength, in Å: the source simulates the band of
        :func:`.sources.wavelength_band` around it, as in ``TRex``. If ``None``, the widest band,
        ``L_MIN`` to ``L_MAX``.
    """

    lambda0: float | None = param(None, "Å", "central wavelength, None for the widest band")


def _source(lambda0: float | None) -> ESS_butterfly:
    """Return the T-Rex source with the band around ``lambda0``, or the widest band if ``None``."""
    return source() if lambda0 is None else source(*wavelength_band(lambda0))


class TRexSource(Instrument):
    """The source of T-Rex and a time and wavelength monitor at the entrance of guide 1.

    The monitor is the entrance of guide 1: its opening, where guide 1 begins, 2.0 mm off the axis
    of the source and turned with the guide. Its window is the widest band, whatever the band of
    the source.

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
        self.source = _source(self.params.lambda0)
        self.source_monitor = source_monitor()

    def update_params(self, **values: Any) -> None:
        """Change control parameters.

        A new ``lambda0`` also creates the source again, with its band.

        Parameters
        ----------
        **values:
            New values by parameter name, e.g. ``update_params(lambda0=2.0)``.

        Raises
        ------
        TypeError
            If a name is not a field of :class:`TRexSourceParams`.
        """
        super().update_params(**values)
        if "lambda0" in values:
            self.source = _source(self.params.lambda0)

    @property
    def instr_params(self) -> list[InstrParam]:
        """McStas instrument parameters: there are none."""
        return []
        # --8<-- [end:trex_source]


def main(
    out: Path | None = None,
    ncount: float = 1e6,
    neutrons: int = 1_000_000,
    seed: int = 1,
    name: str = "trex_source",
    lambda0: float | None = None,
) -> Path:
    """Run both models of :class:`TRexSource`, print their means and save the figure.

    Parameters
    ----------
    out:
        Folder of the figure, ``<name>.png``. If ``None``, the folder of the McStas run, next to
        its ``mccode.h5``.
    ncount:
        Number of McStas neutrons.
    neutrons:
        Number of tof neutrons.
    seed:
        Seed of the random numbers of both models.
    name:
        File name of the figure, without ``.png``.
    lambda0:
        Central wavelength of the band of the source, in Å. If ``None``, the widest band.

    Returns
    -------
    :
        The path of the figure.
    """
    instrument = TRexSource(TRexSourceParams(lambda0=lambda0))

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
    print_means(comparison)
    print(
        f"McStas total {comparison.mcstas_intensity.sum():.4g} n/s, "
        f"{comparison.scale:.4g} n/s per tof neutron"
    )

    fig = plot_comparison(comparison, "T-Rex flux at the entrance of guide 1")
    return save_figure(fig, output_folder(out, folder) / f"{name}.png")


if __name__ == "__main__":
    out_folder = parse_out(__doc__)
    main(out=out_folder, ncount=1e6, neutrons=1_000_000, name="trex_source_1e6")
    # a tenth of the neutrons in both models, for the scale with fewer tof neutrons
    main(out=out_folder, ncount=1e5, neutrons=100_000, name="trex_source_1e5")
    # the band of TRex at its default lambda0, 2 Å, in place of the widest band
    main(out=out_folder, lambda0=2.0, name="trex_source_lambda0_2")
