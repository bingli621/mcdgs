"""The beam at the entrance of the guide system of CSPEC, from McStas and from tof.

:class:`CSPECSource` is the source and a time and wavelength monitor at the entrance of
W03-01-011, the first optics of CSPEC, 1.90424 m from the source and 16.3 mm off its axis, with
the opening of that guide. Nothing stands between the two, so McStas and tof should see the same
pulse there: only their models of the ESS moderator differ.

The monitor is the source monitor of :class:`.CSPEC`, from :func:`.monitors.source_monitor`.

Counts and intensity, scale and error bars are as at the source of T-Rex: McStas sets the
weights of its rays itself, in n/s, while each tof neutron has weight 1, so tof takes its weight
from McStas, the total of McStas over that of tof (:func:`.rescale_tof_to_mcstas`).

:func:`main` runs both models and reads the intensity of each on the bins of the monitor.
:func:`.rescale_tof_to_mcstas` gives the scale of tof, :func:`.mcstas_vs_tof` puts both maps side by
side in n/s, :func:`.print_means` prints their mean time and wavelength, and
:func:`.plot_comparison` draws them. Run the module to do all of it, at 1e6 and at 1e5 neutrons
in both models, and with the band of ``lambda0`` = 3 Å::

    pixi run python -m mcdgs.diagnostics.cspec.source
"""

from dataclasses import dataclass
from pathlib import Path

from mcdgs.core import Instrument, param
from mcdgs.core.components import Arm, ESS_butterfly, Monitor_nD
from mcdgs.diagnostics import cli, compare
from mcdgs.instruments.cspec.guides.frame import arm_iscs
from mcdgs.instruments.cspec.monitors import source_monitor
from mcdgs.instruments.cspec.sources import source, wavelength_band
from mcdgs.mcstas import Simulation, read_mcstas_maps
from mcdgs.reduction.binning import monitor_bins
from mcdgs.tof import build_tof, read_tof_maps


# --8<-- [start:instrument]
@dataclass(frozen=True)
class CSPECSourceParams:
    """Control parameters of :class:`CSPECSource`.

    Parameters
    ----------
    lambda0:
        Central wavelength, in Å: the source simulates the band of
        :func:`.sources.wavelength_band` around it, as in ``CSPEC``. If ``None``, the widest
        band, ``L_MIN`` to ``L_MAX``.
    """

    lambda0: float | None = param(None, "Å", "central wavelength, None for the widest band")


class CSPECSource(Instrument):
    """The source of CSPEC and a time and wavelength monitor at the entrance of W03-01-011.

    The monitor is the entrance of W03-01-011: its opening, where the guide system begins, on
    ``arm_iscs``, the frame of the sheet. Its window is the widest band, whatever the band of
    the source.

    Parameters
    ----------
    params:
        Control parameters. If ``None``, the defaults of :class:`CSPECSourceParams`.
    """

    params_class = CSPECSourceParams
    name = "cspec_source"  # the McStas instrument name

    source: ESS_butterfly
    arm_iscs: Arm
    source_monitor: Monitor_nD

    def __init__(self, params: CSPECSourceParams | None = None) -> None:
        super().__init__(params)
        lambda0 = self.params.lambda0
        self.source = source() if lambda0 is None else source(*wavelength_band(lambda0))
        self.arm_iscs = arm_iscs()  # the frame the entrance of the guide system is placed in
        self.source_monitor = source_monitor()


def define_instrument(lambda0: float | None = None) -> CSPECSource:
    """Return the instrument of this diagnostic.

    Parameters
    ----------
    lambda0:
        Central wavelength of the band of the source, in Å. If ``None``, the widest band, 0.1
        to 10 Å.

    Returns
    -------
    :
        :class:`CSPECSource` with that band.
    """
    source_params = CSPECSourceParams(lambda0=lambda0)
    return CSPECSource(source_params)


# --8<-- [end:instrument]


def main(
    out: Path | None = None,
    ncount: float = 1e6,
    neutrons: int = 1_000_000,
    seed: int = 1,
    lambda0: float | None = None,
    folder_name: str = "cspec_source",
) -> Path:
    """Run both models of :class:`CSPECSource`, print their means and save the figure.

    Parameters
    ----------
    out:
        Folder of the figure, ``<folder_name>.png``. If ``None``, ``cspec_source`` in the
        results folder of :meth:`.Workspace.find`.
    ncount:
        Number of McStas neutrons.
    neutrons:
        Number of tof neutrons.
    seed:
        Seed of the random numbers of both models.
    lambda0:
        Central wavelength of the band of the source, in Å. If ``None``, the widest band.
    folder_name:
        Name of the McStas run folder, and file name of the figure, without ``.png``.

    Returns
    -------
    :
        The path of the figure.
    """
    instrument = define_instrument(lambda0)
    # --8<-- [start:run_mcstas]
    # run McStas
    sim = Simulation(instrument)
    mcstas_run = sim.run(folder_name=folder_name, ncount=ncount, seed=seed)
    # --8<-- [end:run_mcstas]
    # --8<-- [start:run_tof]
    # run tof, then read the monitor of both on the bins of the McStas map: time, wavelength
    tof_model = build_tof(instrument, neutrons=neutrons, seed=seed)
    tof_run = tof_model.run()
    # --8<-- [end:run_tof]

    # the bins of the map of the monitor, time and wavelength, read from the monitor
    bins = {"source_monitor": monitor_bins(instrument.source_monitor.bin_limits)}
    # the map of each model on those bins, a scipp data array with its errors as variances:
    # McStas in n/s, tof in neutrons, a bin of N neutrons with the Poisson error √N
    mcstas_map = read_mcstas_maps(mcstas_run, bins)["source_monitor"]
    tof_counts = read_tof_maps(tof_run, bins)["source_monitor"]

    # tof has no flux: each of its neutrons stands for the total of McStas over that of tof
    scale = compare.rescale_tof_to_mcstas(tof_counts, mcstas_map)  # n/s per tof neutron
    # both maps side by side in n/s: tof times the scale, its Poisson errors scaled along
    comparison = compare.mcstas_vs_tof(
        "source_monitor", mcstas_map, tof_counts, ncount, neutrons, scale
    )

    # the mean time and wavelength of each model, weighted by its intensity, and the scale
    compare.print_means(comparison)
    print(f"McStas total {mcstas_map.sum().value:.4g} n/s, {scale:.4g} n/s per tof neutron")

    # top row: the time and wavelength map of each model, on one color scale; bottom row: the
    # wavelength spectrum and the time profile of both, with their error bars
    fig = compare.plot_comparison(comparison, "CSPEC flux at the entrance of the guide system")

    return cli.save_figure(fig, cli.output_folder(out, "cspec_source") / f"{mcstas_run.name}.png")


if __name__ == "__main__":
    out_folder = cli.parse_out(__doc__)
    main(out=out_folder, ncount=1e6, neutrons=1_000_000, folder_name="cspec_source_1e6")
    # a tenth of the neutrons in both models, for the scale with fewer tof neutrons
    main(out=out_folder, ncount=1e5, neutrons=100_000, folder_name="cspec_source_1e5")
    # the band of CSPEC at its default lambda0, 3 Å, in place of the widest band
    main(out=out_folder, lambda0=3.0, folder_name="cspec_source_lambda0_3A")
