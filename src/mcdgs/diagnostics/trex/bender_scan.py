"""How the bender passes each band of wavelengths, from a scan of its rotation.

:class:`TRexBenderScan` is T-Rex up to the beam monitor after guide 4, with the bender in, as the
:class:`.Extraction` has it. The scan reads its beam monitor: its time and wavelength map, summed
over time, in bands of :data:`BAND_WIDTH` centered on its multiples.
:func:`run_bender_scan` scans ``b_rot`` with the bender in, and runs once with it out;
:func:`plot_bender_scan` reads the beam monitor of both, the scan with :func:`.read_scan_maps` and
the run with ``mcstasscript.load_data`` and :func:`.read_mcstas_map`, sums each map over time and into the bands, and draws the
intensity per band of the scan over that of the same band with the bender out. :func:`main`
does both; run the module for it::

    pixi run python -m mcdgs.diagnostics.trex.bender_scan
"""

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import mcstasscript as ms
import numpy as np
import scipp as sc
from matplotlib.figure import Figure
from numpy.typing import NDArray

from mcdgs.core import InstrParam, Instrument, param
from mcdgs.core.components import (
    Arm,
    ESS_butterfly,
    Guide_gravity,
    Monitor_nD,
    Pol_bender,
    Progress_bar,
)
from mcdgs.diagnostics import cli
from mcdgs.instruments.trex.guides.bispectral import (
    BENDER_ROTATION,
    arm_guide_4,
    bender,
    guide_4a,
    guide_4b,
)
from mcdgs.instruments.trex.guides.extraction import guides_extraction
from mcdgs.instruments.trex.monitors import beam_monitor, source_monitor
from mcdgs.instruments.trex.sources import source
from mcdgs.mcstas import Run, Simulation, Workspace, read_mcstas_map, read_scan_maps
from mcdgs.reduction.binning import monitor_bins

BAND_WIDTH = 0.5  # Å: the bands are centered on its multiples, 0.5, 1.0, 1.5 Å, ...
BAND_LIMITS = (0.25, 9.75)  # Å: the edges of the first and the last band, 0.5 and 9.5 Å


# --8<-- [start:instrument]
@dataclass(frozen=True)
class TRexBenderScanParams:
    """Control parameters of :class:`TRexBenderScan`.

    Parameters
    ----------
    bender_in:
        True uses the bender and the cold guide 4b, False uses the thermal guide 4a.
    b_rot:
        Rotation of the bender about +y, in deg.
    """

    bender_in: bool = param(True, "", "True = bender in, False = bender out")
    b_rot: float = param(BENDER_ROTATION, "deg", "rotation of the cold neutron bender")


class TRexBenderScan(Instrument):
    """T-Rex up to the beam monitor after guide 4, with the bender in by default.

    The source keeps the widest band, 0.1 to 10 Å. The beam monitor covers the bands,
    :data:`BAND_LIMITS`, so that the edge of every band is also the edge of a bin. Guide 4a, the
    bender and guide 4b are all in the instrument; ``bender_in`` chooses the path, as an
    instrument parameter, so that a run of the scan can set it.

    Parameters
    ----------
    params:
        Control parameters. If ``None``, the defaults of :class:`TRexBenderScanParams`.
    """

    params_class = TRexBenderScanParams
    name = "trex_bender_scan"  # the McStas instrument name

    origin: Progress_bar
    source: ESS_butterfly
    source_monitor: Monitor_nD
    arm_nboa: Arm
    guide_1: Guide_gravity
    guide_2: Guide_gravity
    guide_3: Guide_gravity
    arm_guide_4: Arm
    guide_4a: Guide_gravity
    bender: Pol_bender
    guide_4b: Guide_gravity
    beam_monitor: Monitor_nD

    def __init__(self, params: TRexBenderScanParams | None = None) -> None:
        super().__init__(params)
        self.origin = Progress_bar(at=[0, 0, 0])
        self.source = source()  # the widest band
        self.source_monitor = source_monitor()
        self.add_components(guides_extraction())
        self.arm_guide_4 = arm_guide_4()
        self.guide_4a = guide_4a(when="bender_in==0")
        self.bender = bender(when="bender_in==1")
        self.guide_4b = guide_4b(when="bender_in==1")
        self.beam_monitor = beam_monitor(lambda_limits=BAND_LIMITS)  # on the edges of the bands

    @property
    def instr_params(self) -> list[InstrParam]:
        """The McStas instrument parameters: ``bender_in`` as 0 or 1, and ``b_rot``."""
        return [
            InstrParam("bender_in", int(self.params.bender_in), "", "bender in = 1, out = 0"),
            InstrParam("b_rot", self.params.b_rot, "deg", "rotation of the cold neutron bender"),
        ]


def define_instrument() -> TRexBenderScan:
    """Return the instrument of this diagnostic: :class:`TRexBenderScan`, the bender in.

    The runs of the scan still set ``bender_in`` themselves: the scan with the bender in, the run
    without it with it out.
    """
    params = TRexBenderScanParams()  # the default control parameters: the bender in
    return TRexBenderScan(params)


# --8<-- [end:instrument]


def run_bender_scan(
    instrument: TRexBenderScan,
    angles: NDArray | None = None,
    ncount: float = 1e6,
    seed: int = 1,
    workspace: Workspace | None = None,
) -> tuple[Run, Run]:
    """Scan the rotation of the bender, and run once without it.

    The scan goes to the run folder ``bender_scan`` of the workspace, the run without the
    bender to ``bender_out``; each replaces the one of an earlier call.

    Parameters
    ----------
    instrument:
        The instrument to run, as :func:`define_instrument` gives it. The runs set
        ``bender_in`` themselves.
    angles:
        The values of ``b_rot``, in deg. If ``None``, 19 from -0.8 to 0.1 deg.
    ncount:
        Number of McStas neutrons, per angle and in the run without the bender.
    seed:
        Seed of the random numbers.
    workspace:
        The folders to write into. If ``None``, :meth:`.Workspace.find`.

    Returns
    -------
    :
        The scan and the run without the bender.
    """
    if angles is None:
        angles = np.linspace(-0.8, 0.1, 19)

    # --8<-- [start:run_scan]
    sim = Simulation(instrument, workspace)
    bender_in_scan = sim.scan(
        variable="b_rot",
        values=angles,
        folder_name="bender_scan",
        ncount=ncount,
        seed=seed,
        bender_in=1,
    )
    print(f"Saved {bender_in_scan.data_file}")
    bender_out_run = sim.run(
        folder_name="bender_out",
        ncount=ncount,
        seed=seed,
        bender_in=0,
    )
    # --8<-- [end:run_scan]
    print(f"Saved {bender_out_run.data_file}")
    return bender_in_scan, bender_out_run


def plot_bender_scan(
    instrument: TRexBenderScan, bender_in_scan: Run, bender_out_run: Run
) -> Figure:
    """Draw the scan of :func:`run_bender_scan`.

    Each curve is one wavelength band, divided by the same band without the bender, so the
    bender-out intensity is 1. The bands take matplotlib's default colors in turn, solid for
    the first ten and dash-dotted for the next, and the legend gives the center and half width
    of each.

    Parameters
    ----------
    instrument:
        The instrument of the runs, whose beam monitor gives the bins of its maps.
    bender_in_scan:
        The scan with the bender in.
    bender_out_run:
        The run with the bender out.

    Returns
    -------
    :
        The figure.
    """
    bins = monitor_bins(instrument.beam_monitor.bin_limits)
    # --8<-- [start:read]
    bender_in = read_scan_maps("beam_monitor", bender_in_scan, bins)
    bender_out = read_mcstas_map("beam_monitor", ms.load_data(bender_out_run.folder), bins)
    # over time, then over the bins of each band: the edges of the bands, 0.25, 0.75, ...
    # 9.75 Å, are edges of bins, so each bin goes whole to one band, and the variances add
    bands = sc.linspace("wavelength", *BAND_LIMITS, 20, unit="Å")
    bender_in = bender_in.sum("toa").rebin(wavelength=bands)
    bender_out = bender_out.sum("toa").rebin(wavelength=bands)
    # --8<-- [end:read]
    # --8<-- [start:transmission]
    # scipp does not broadcast variances, which would correlate the points: each point gets a
    # copy of the reference of its own, and scipp adds the relative errors in quadrature
    ratio = bender_in / bender_out.broadcast(sizes=bender_in.sizes).copy()
    # --8<-- [end:transmission]
    angles = ratio.coords["b_rot"].values
    centers = sc.midpoints(ratio.coords["wavelength"]).values

    fig, ax = plt.subplots(figsize=(8, 4.8), layout="constrained")
    n_colors = 10  # the default color cycle, C0 to C9: past it the colors come again
    for band, center in enumerate(centers):
        curve = ratio["wavelength", band]
        ax.errorbar(
            angles,
            curve.values,
            sc.stddevs(curve).values,
            color=f"C{band % n_colors}",
            linestyle="-" if band < n_colors else "-.",
            label=f"{center:.2f}±{BAND_WIDTH / 2:.2f} Å",
        )
    ax.axhline(1, color="black", linestyle="--", label="bender out")
    ax.axvline(
        BENDER_ROTATION, color="black", linestyle=":", label=f"b_rot = {BENDER_ROTATION:.3f} deg"
    )
    ax.set_ylim(bottom=0)
    ax.set_xlabel("b_rot [deg]")
    ax.set_ylabel("Intensity at the beam monitor / bender out")
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize="small")  # beside the plot
    return fig


def main(out: Path | None = None, ncount: float = 1e6, seed: int = 1) -> Path:
    """Run the scan of the bender and save its figure.

    Parameters
    ----------
    out:
        Folder of the figures. If ``None``, ``trex_bender_scan`` in the results folder of
        :meth:`.Workspace.find`.
    ncount:
        Number of McStas neutrons, per angle and in the run without the bender.
    seed:
        Seed of the random numbers.

    Returns
    -------
    :
        The path of the figure, ``bender_scan.png``.
    """
    instrument = define_instrument()
    bender_in_scan, bender_out_run = run_bender_scan(instrument, ncount=ncount, seed=seed)
    fig = plot_bender_scan(instrument, bender_in_scan, bender_out_run)
    return cli.save_figure(fig, cli.output_folder(out, "trex_bender_scan") / "bender_scan.png")


if __name__ == "__main__":
    main(out=cli.parse_out(__doc__))
