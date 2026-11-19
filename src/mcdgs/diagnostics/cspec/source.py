"""The beam at the entrance of the guide system of CSPEC, from McStas and from tof.

:class:`CSPECSource` is the source and a time and wavelength monitor at the entrance of
W03-01-011, the first optics of CSPEC, 1.90424 m from the source and 16.3 mm off its axis, with
the opening of that guide. Nothing stands between the two, so McStas and tof should see the same
pulse there: only their models of the ESS moderator differ.

The model of CSPEC has no monitors: :func:`source_monitor` is this diagnostic's own, placed as
the entrance of W03-01-011 is in :class:`.CSPEC`.

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

import mcstasscript as ms

from mcdgs.core import Instrument, param
from mcdgs.core.components import Arm, ESS_butterfly, Monitor_nD
from mcdgs.core.components.monitors import time_lambda_options
from mcdgs.diagnostics import cli, compare
from mcdgs.instruments.cspec.guides.bender import (
    GUIDE_01_011_IN_HEIGHT,
    GUIDE_01_011_IN_WIDTH,
    GUIDE_START_Z,
    guides_bender,
)
from mcdgs.instruments.cspec.guides.frame import MM_PER_M, arm_iscs
from mcdgs.instruments.cspec.sources import source, wavelength_band
from mcdgs.mcstas import Simulation, read_mcstas_map
from mcdgs.reduction.binning import monitor_bins
from mcdgs.tof import build_tof, read_tof_map

SOURCE_MONITOR_LAMBDA_RESOLUTION = 0.05  # Å, the width of a wavelength bin, as for T-Rex
SOURCE_MONITOR_TIME_RESOLUTION = 5e-5  # s, the width of a time bin, 0.05 ms, as for T-Rex


def source_monitor(
    lambda_resolution: float = SOURCE_MONITOR_LAMBDA_RESOLUTION,
    time_resolution: float = SOURCE_MONITOR_TIME_RESOLUTION,
) -> Monitor_nD:
    """Create the time and wavelength monitor at the entrance of W03-01-011.

    It is the entrance of W03-01-011: its opening, at the place and with the rotation of that
    guide on ``arm_iscs``, so it counts what enters the guide system. It records one 2D map of
    time against wavelength, and leaves the neutrons unchanged. Its axes are those of
    :func:`.time_lambda_options` at ``GUIDE_START_Z``, without limits: the widest band the
    source simulates, ``L_MIN`` to ``L_MAX``, and the time window in which it reaches the
    monitor.

    Parameters
    ----------
    lambda_resolution:
        Width of a wavelength bin, in Å. The axis has the fewest bins that cover the band, see
        :func:`.bins_from_limits`.
    time_resolution:
        Width of a time bin, in s, likewise.

    Returns
    -------
    :
        A new monitor, relative to ``arm_iscs``.
    """
    to_bwc1, _, _ = guides_bender()
    entrance = to_bwc1["guide_01_011"]
    return Monitor_nD(
        xwidth=GUIDE_01_011_IN_WIDTH / MM_PER_M,
        yheight=GUIDE_01_011_IN_HEIGHT / MM_PER_M,
        restore_neutron=1,
        options=time_lambda_options(GUIDE_START_Z, None, None, lambda_resolution, time_resolution),
        filename="source_monitor",  # the name of the monitor, in both models
        at=list(entrance.at),
        rotated=list(entrance.rotated),
        relative=entrance.relative,
    )


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
    # run McStas and load its monitors, in n/s per bin
    sim = Simulation(instrument)
    run = sim.run(folder_name=folder_name, ncount=ncount, seed=seed)
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
    scale = compare.rescale_tof_to_mcstas(tof_counts, mcstas_map)  # n/s per tof neutron
    # both maps side by side in n/s: tof times the scale, its Poisson errors scaled along
    comparison = compare.mcstas_vs_tof(
        "source_monitor", mcstas_map, tof_counts, ncount, neutrons, scale
    )
    # --8<-- [end:scale_tof]

    # the mean time and wavelength of each model, weighted by its intensity, and the scale
    compare.print_means(comparison)
    print(f"McStas total {mcstas_map.sum().value:.4g} n/s, {scale:.4g} n/s per tof neutron")

    # --8<-- [start:plot]
    # top row: the time and wavelength map of each model, on one color scale; bottom row: the
    # wavelength spectrum and the time profile of both, with their error bars
    fig = compare.plot_comparison(comparison, "CSPEC flux at the entrance of the guide system")
    # --8<-- [end:plot]

    return cli.save_figure(fig, cli.output_folder(out, "cspec_source") / f"{run.name}.png")


if __name__ == "__main__":
    out_folder = cli.parse_out(__doc__)
    main(out=out_folder, ncount=1e6, neutrons=1_000_000, folder_name="cspec_source_1e6")
    # a tenth of the neutrons in both models, for the scale with fewer tof neutrons
    main(out=out_folder, ncount=1e5, neutrons=100_000, folder_name="cspec_source_1e5")
    # the band of CSPEC at its default lambda0, 3 Å, in place of the widest band
    main(out=out_folder, lambda0=3.0, folder_name="cspec_source_lambda0_3A")
