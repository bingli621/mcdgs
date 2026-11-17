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

    pixi run python -m mcdgs.diagnostics.trex.source
"""

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import mcstasscript as ms

from mcdgs.core import Instrument, param
from mcdgs.core.components import ESS_butterfly, Monitor_nD
from mcdgs.diagnostics.cli import output_folder, parse_out, save_figure
from mcdgs.diagnostics.compare import intensity_label, mean_coordinate, plot_map, plot_profile
from mcdgs.instruments.trex.monitors import source_monitor
from mcdgs.instruments.trex.sources import source, wavelength_band
from mcdgs.mcstas import Simulation, read_mcstas_map
from mcdgs.reduction.binning import monitor_bins
from mcdgs.tof import build_tof, read_tof_map


# --8<-- [start:instrument]
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
        lambda0 = self.params.lambda0
        self.source = source() if lambda0 is None else source(*wavelength_band(lambda0))
        self.source_monitor = source_monitor()


def define_instrument(lambda0: float | None = None) -> TRexSource:
    """Return the instrument of this diagnostic.

    Parameters
    ----------
    lambda0:
        Central wavelength of the band of the source, in Å. If ``None``, the widest band, 0.1
        to 10 Å.

    Returns
    -------
    :
        :class:`TRexSource` with that band.
    """
    source_params = TRexSourceParams(lambda0=lambda0)
    return TRexSource(source_params)


# --8<-- [end:instrument]


def main(
    out: Path | None = None,
    ncount: float = 1e6,
    neutrons: int = 1_000_000,
    seed: int = 1,
    lambda0: float | None = None,
    folder_name: str = "trex_source",
) -> Path:
    """Run both models of :class:`TRexSource`, print their means and save the figure.

    Parameters
    ----------
    out:
        Folder of the figure, ``<name>.png``. If ``None``, ``trex_source`` in the results folder of
        :meth:`.Workspace.find`.
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
    instrument = define_instrument(lambda0)
    # --8<-- [start:run_mcstas]
    # run McStas and load its monitors, in n/s per bin
    sim = Simulation(instrument)
    run = sim.run(folder_name=folder_name, ncount=ncount)
    mcstas_data = ms.load_data(run.folder)
    # --8<-- [end:run_mcstas]
    # --8<-- [start:run_tof]
    # run tof, then read the monitor of both on the bins of the McStas map: time, wavelength
    tof_model = build_tof(instrument, neutrons=neutrons, seed=seed)
    tof_result = tof_model.run()
    # --8<-- [end:run_tof]

    # --8<-- [start:read_maps]
    # the bins of the map of the monitor, time and wavelength, read from the monitor
    bins = monitor_bins(instrument.source_monitor.bin_limits)
    # the map of each model on those bins, a scipp data array with its errors as variances:
    # McStas in n/s, tof in neutrons, a bin of N neutrons with the Poisson error √N
    mcstas_map = read_mcstas_map("source_monitor", mcstas_data, bins)
    tof_counts = read_tof_map("source_monitor", tof_result, bins)
    # --8<-- [end:read_maps]

    # --8<-- [start:scale_tof]
    # tof has no flux: each of its neutrons stands for the total of McStas over that of tof
    scale = mcstas_map.sum().value / tof_counts.sum().value  # n/s per tof neutron
    tof_map = tof_counts * scale  # the errors are scaled along
    tof_map.unit = "1/s"
    # --8<-- [end:scale_tof]

    # the mean time and wavelength of each model, weighted by its intensity
    for label, data in (("McStas", mcstas_map), ("tof", tof_map)):
        mean_time = mean_coordinate(data, "toa").to(unit="ms").value
        mean_wavelength = mean_coordinate(data, "wavelength").to(unit="Å").value
        print(f"{label:>6}: mean time {mean_time:.3f} ms, mean wavelength {mean_wavelength:.3f} Å")
    print(f"McStas total {mcstas_map.sum().value:.4g} n/s, {scale:.4g} n/s per tof neutron")

    # --8<-- [start:plot]
    fig, ((mcstas_ax, tof_ax), (spectrum_ax, profile_ax)) = plt.subplots(
        2, 2, figsize=(12, 9), layout="constrained"
    )
    fig.suptitle("T-Rex flux at the entrance of guide 1")

    # top row: the time and wavelength map of each model, on one color scale
    vmax = max(mcstas_map.max().value, tof_map.max().value)
    plot_map(mcstas_ax, mcstas_map, vmax=vmax)
    mesh = plot_map(tof_ax, tof_map, vmax=vmax)
    mcstas_ax.set_title("McStas")
    tof_ax.set_title("tof")
    fig.colorbar(mesh, ax=[mcstas_ax, tof_ax], label=intensity_label(mcstas_map))

    # bottom row: the wavelength spectrum, each map summed over time, and the time profile,
    # summed over wavelength. Summing adds the variances, so the errors add in quadrature
    mcstas_label = f"McStas, ncount = {int(ncount):,}"
    tof_label = f"tof, {neutrons:,} neutrons × {scale:.3g} n/s each"
    for ax, summed in ((spectrum_ax, "toa"), (profile_ax, "wavelength")):
        plot_profile(ax, mcstas_map.sum(summed), mcstas_label)
        plot_profile(ax, tof_map.sum(summed), tof_label)
        ax.set_ylim(bottom=0)
        ax.legend()
    # --8<-- [end:plot]

    return save_figure(fig, output_folder(out, "trex_source") / f"{run.name}.png")


if __name__ == "__main__":
    out_folder = parse_out(__doc__)
    main(out=out_folder, ncount=1e6, neutrons=1_000_000, folder_name="trex_source_1e6")
    # a tenth of the neutrons in both models, for the scale with fewer tof neutrons
    main(out=out_folder, ncount=1e5, neutrons=100_000, folder_name="trex_source_1e5")
    # the band of TRex at its default lambda0, 2 Å, in place of the widest band
    main(out=out_folder, lambda0=2.0, folder_name="trex_source_lambda0_2A")
