"""T-Rex up to guide 4, the thermal and the cold path, to study the extraction guides.

:class:`Extraction` is the section from the source to the beam monitor after guide 4, at 5.9 m,
with the same components as :class:`.TRex`, but for the band: its source and its beam monitor
keep the widest band, 0.1 to 10 Å, where those of :class:`.TRex` follow ``lambda0``. It also has
the source monitor of :class:`.TRexSource`, at the entrance of guide 1.

:func:`diagnose_path` compares McStas and tof at the beam monitor, on one path through guide 4.
tof has no guides: every neutron of its pulse reaches the monitor. tof is scaled to the total of
McStas at the source monitor, where both see the same beam, and by ``transmission`` on top of
it, so the beam monitor shows what the guides let through in McStas against that factor. The
McStas runs also draw McStasScript's diagnostics, the beam after the source, after guide 3 and
at the exit of guide 4b, and its diagram of the intensity before each component. :func:`main`
draws the extraction from above and runs the path without the bender, through guide 4a; the
path with it is run by :func:`.bender_scan.main`. Run the module for it::

    pixi run python -m mcdgs.diagnostics.trex.extraction
"""

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from unittest.mock import patch

import matplotlib.pyplot as plt
import mcstasscript as ms
from matplotlib.figure import Figure

from mcdgs.core import InstrParam, Instrument, param
from mcdgs.core.components import (
    Arm,
    ESS_butterfly,
    Guide_gravity,
    Monitor_nD,
    Pol_bender,
    Progress_bar,
)
from mcdgs.diagnostics.cli import output_folder, parse_out, save_figure
from mcdgs.diagnostics.compare import plot_comparison, read_scaled_comparisons
from mcdgs.diagnostics.drawing import draw_instrument
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
from mcdgs.mcstas import Run, RunType, Workspace, build_mcstas, recorded_instr_params
from mcdgs.reduction.binning import monitor_bins
from mcdgs.tof import build_tof

# McStas neutrons of the run with McStasScript's diagnostics, which lists every ray at each of its
# three positions: about 180 MB at 1e6
DIAGNOSTICS_NCOUNT = 1e6

# n/s: the range of the intensity axis of McStasScript's diagram, from the source to the beam
# monitor on both paths, about 1e11 to 1e12 n/s at the defaults
DIAGRAM_INTENSITY_LIMITS = (5e10, 1.5e12)

# m: the view from above, z along the beam from the source to the entrance of guide 5 (cell
# B32), just after the beam monitor, then x across it, wide enough for the moderator
TOP_VIEW_LIMITS = (0.0, 5.916, -0.15, 0.15)


# --8<-- [start:extraction]
@dataclass(frozen=True)
class ExtractionParams:
    """Control parameters of the extraction section.

    Frozen: change them with ``Extraction.update_params``. Each field's unit and comment are in
    its metadata.

    Parameters
    ----------
    bender_in:
        True uses the bender and the cold guide 4b, False uses the thermal guide 4a.
    b_rot:
        Rotation of the bender about +y, in deg.
    """

    bender_in: bool = param(False, "", "True = bender in, False = bender out")
    b_rot: float = param(BENDER_ROTATION, "deg", "rotation of the cold neutron bender")


DEFAULT_PARAMS = ExtractionParams()  # frozen, so one shared default is safe


class Extraction(Instrument):
    """T-Rex up to guide 4, the thermal and the cold path, and the beam monitor after it.

    Parameters
    ----------
    params:
        Control parameters.
    """

    params_class = ExtractionParams
    name = "trex_extraction"  # the McStas instrument name

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

    def __init__(self, params: ExtractionParams = DEFAULT_PARAMS) -> None:
        super().__init__(params)
        self.origin = Progress_bar(at=[0, 0, 0])
        self.source = source()  # the widest band, as the beam monitor's default limits
        self.source_monitor = source_monitor()
        self.add_components(guides_extraction())
        self.arm_guide_4 = arm_guide_4()
        self.guide_4a = guide_4a(when="bender_in==0")
        self.bender = bender(when="bender_in==1")
        self.guide_4b = guide_4b(when="bender_in==1")
        self.beam_monitor = beam_monitor()

    @property
    def instr_params(self) -> list[InstrParam]:
        """The McStas instrument parameters: ``bender_in`` as 0 or 1, and ``b_rot``."""
        return [
            InstrParam("bender_in", int(self.params.bender_in), "", "bender in = 1, out = 0"),
            InstrParam("b_rot", self.params.b_rot, "deg", "rotation of the cold neutron bender"),
        ]
        # --8<-- [end:extraction]


def _zoom_intensity(diagram: Figure, limits: tuple[float, float]) -> None:
    """Zoom the intensity axis of McStasScript's diagram, and its ray count with it.

    ``show_diagram`` takes no range for the intensity: its ``limits`` are those of a
    ``variable``. The ray count keeps its ratio to the intensity, so both share their ticks.
    """
    intensity_ax = next(ax for ax in diagram.axes if ax.get_xlabel() == "Intensity [n/s]")
    rays_ax = next(ax for ax in diagram.axes if ax.get_xlabel() == "Ray count")
    rays_per_intensity = rays_ax.get_xlim()[1] / intensity_ax.get_xlim()[1]
    intensity_ax.set_xlim(limits)
    rays_ax.set_xlim(limits[0] * rays_per_intensity, limits[1] * rays_per_intensity)


def save_run_record(
    folder: Path,
    extraction: Extraction,
    started: datetime,
    bender_in: bool,
    ncount: float,
    seed: int,
) -> Run:
    """Save the record of a run of the extraction that McStasScript started itself.

    The runs of :func:`diagnose_path` go through McStasScript's diagnostics and its diagram,
    not through a :class:`.Simulation`, which would save the record. With it the run is found
    and cleaned as any other. McStasScript writes them in McStas's own format, so the
    ``data_file`` of the record, a ``mccode.h5``, does not exist: the folder is loaded with
    ``ms.load_data``.

    Parameters
    ----------
    folder:
        The run folder, written by the run.
    extraction:
        The instrument.
    started:
        When the run started.
    bender_in:
        The path of the run: True for the bender and guide 4b.
    ncount, seed:
        The number of neutrons and the seed of the run.

    Returns
    -------
    :
        The run, finished now and saved in its folder.
    """
    run = Run(
        folder=folder,
        instrument=extraction.name,
        run_type=RunType.RUN,
        started=started,
        instr_params=recorded_instr_params(extraction.instr_params, {"bender_in": int(bender_in)}),
        ncount=ncount,
        seed=seed,
    )
    return run.finish()


def diagnose_path(
    bender_in: bool,
    out: Path | None = None,
    ncount: float = 1e7,
    neutrons: int = 1_000_000,
    seed: int = 1,
    transmission: float = 0.2,
    diagnostics_ncount: float = DIAGNOSTICS_NCOUNT,
) -> list[Path]:
    """Run one path through guide 4 in McStas and tof, and save its three figures.

    McStas runs twice, once with McStasScript's diagnostics, the beam after the source, after
    guide 3 and at the exit of guide 4b, and once for its diagram of the intensity before each
    component, which also gives the monitors; tof runs once. The diagnostics list every ray at
    their three positions, so their run has fewer neutrons. tof is scaled to the total of
    McStas at the source monitor, where both see the same beam, and by ``transmission`` on top
    of it at the beam monitor. The figures, the diagnostics, the diagram and the comparison, are
    ``trex_extraction_bender_out_diagnose.png``, ``trex_extraction_bender_out_intensity.png``
    and ``trex_extraction_bender_out.png``, or ``..._in``.

    Parameters
    ----------
    bender_in:
        True for the bender and guide 4b, False for guide 4a.
    out:
        Folder of the figures. If ``None``, ``trex_extraction`` in the results folder of
        :meth:`.Workspace.find`.
    ncount:
        Number of McStas neutrons of the run of the diagram and the monitors.
    neutrons:
        Number of tof neutrons.
    seed:
        Seed of the random numbers of both models.
    transmission:
        The factor on tof at the beam monitor, on top of the scale of the source monitor: the
        fraction of the beam the guides would let through, which tof does not model. 0.2 is
        an estimate, close to the 0.17 McStas has with the bender out at the defaults.
    diagnostics_ncount:
        Number of McStas neutrons of the run with the diagnostics.

    Returns
    -------
    :
        The paths of the diagnostics, the diagram and the comparison.
    """
    extraction = Extraction()
    path = "bender in, guide 4b" if bender_in else "bender out, guide 4a"
    run_name = f"trex_extraction_bender_{'in' if bender_in else 'out'}"

    # the bin edges of each monitor, read from it so that they cannot drift apart
    source_bins = monitor_bins(extraction.source_monitor.bin_limits)
    bins = monitor_bins(extraction.beam_monitor.bin_limits)

    # tof has no bender: its run is the same on both paths
    tof_result = build_tof(extraction, neutrons=neutrons, seed=seed).run()

    # --8<-- [start:diagnostics]
    # McStasScript's diagnostics: the beam after the source, guide 3 and guide 4b, its profile,
    # wavelength against x, and its spectrum
    workspace = Workspace.find()
    mcstas_instr = build_mcstas(extraction, input_path=workspace.build_folder(extraction.name))
    diag = ms.Diagnostics(mcstas_instr)  # works on a copy: mcstas_instr is left as it is
    for component in ("source", "guide_3", "guide_4b"):
        # rays="all": a limited list stops early, which leaves the intensities too low
        diag.add_point(after=component, rays="all")  # pyright: ignore[reportArgumentType]
    diag.add_view("x", "y", bins=100, colormap="turbo")  # McStasScript's default is jet
    diag.add_view("x", "l", bins=100, colormap="turbo")
    diag.add_view("l")
    # --8<-- [end:diagnostics]

    started = datetime.now().astimezone()
    folder = workspace.new_run_folder(run_name, "diagnose")
    diag.set_parameters(bender_in=int(bender_in))
    diag.settings(
        ncount=diagnostics_ncount,
        suppress_output=True,
        mpi="auto",
        seed=seed,
        output_path=str(folder),
    )
    diag.run()
    save_run_record(folder, extraction, started, bender_in, diagnostics_ncount, seed)
    diag.plot()  # draws on a new figure, which it does not return
    figure_folder = output_folder(out, "trex_extraction")
    figures = [save_figure(plt.gcf(), figure_folder / f"{run_name}_diagnose.png", dpi=100)]

    # McStasScript's diagram of the components, with the intensity and the rays before each:
    # a McStas run of its own, in a folder of its own, which also holds the monitors
    started = datetime.now().astimezone()
    intensity_folder = workspace.new_run_folder(f"{run_name}_intensity", "diagram")
    mcstas_instr.set_parameters(bender_in=int(bender_in))
    mcstas_instr.settings(
        ncount=ncount,
        suppress_output=True,
        mpi="auto",
        seed=seed,
        output_path=str(intensity_folder),
    )
    # --8<-- [start:diagram]
    # show_diagram ends in plt.show(), which in a window waits for it to close, and the figure
    # goes with it: hold it back, to save the figure as the others
    with patch.object(plt, "show"):
        mcstas_instr.show_diagram(analysis=True)
    diagram = plt.gcf()
    _zoom_intensity(diagram, DIAGRAM_INTENSITY_LIMITS)
    # --8<-- [end:diagram]
    save_run_record(intensity_folder, extraction, started, bender_in, ncount, seed)
    figures.append(save_figure(diagram, figure_folder / f"{run_name}_intensity.png", dpi=100))

    # --8<-- [start:scale_tof]
    # tof's counts are scaled to the total of McStas at the source monitor, in n/s, and by the
    # transmission on top of it at the beam monitor
    mcstas_data = ms.load_data(str(intensity_folder))
    comparisons = read_scaled_comparisons(
        mcstas_data,
        tof_result,
        {"source_monitor": source_bins, "beam_monitor": bins},
        ncount,
        neutrons,
        transmission=transmission,
    )
    at_source, comparison = comparisons["source_monitor"], comparisons["beam_monitor"]
    # --8<-- [end:scale_tof]
    mcstas_total = comparison.mcstas_intensity.sum()
    print(
        f"{path}: McStas total {mcstas_total:.4g} n/s, "
        f"{mcstas_total / at_source.mcstas_intensity.sum():.4g} of that at the source, "
        f"{comparison.scale:.4g} n/s per tof neutron"
    )
    fig = plot_comparison(comparison, f"T-Rex flux at the beam monitor after guide 4, {path}")
    figures.append(save_figure(fig, figure_folder / f"{run_name}.png"))
    return figures


def main(
    out: Path | None = None,
    ncount: float = 1e7,
    neutrons: int = 1_000_000,
    seed: int = 1,
    transmission: float = 0.2,
) -> list[Path]:
    """Draw the extraction from above, and run the path without the bender.

    The view from above is ``trex_extraction_top.png``; the path through guide 4a is
    :func:`diagnose_path` with the bender out. The path with the bender in is drawn by
    :func:`.bender_scan.main`.

    Parameters
    ----------
    out:
        Folder of the figures. If ``None``, ``trex_extraction`` in the results folder of
        :meth:`.Workspace.find`.
    ncount, neutrons, seed, transmission:
        As in :func:`diagnose_path`.

    Returns
    -------
    :
        The paths of the view from above, then of the diagnostics, the diagram and the
        comparison.
    """
    extraction = Extraction()
    build_folder = Workspace.find().build_folder(extraction.name)
    # the geometry: mcdisplay draws every component, whatever its WHEN
    mcstas_instr = build_mcstas(extraction, input_path=build_folder)
    figure_folder = output_folder(out, "trex_extraction")
    # from above: z along the beam, x across it, 16 x 4 in, within TOP_VIEW_LIMITS
    top = draw_instrument(
        mcstas_instr,
        figure_folder / "trex_extraction_top.png",
        equal_aspect=False,
        projection="zx",
        limits=list(TOP_VIEW_LIMITS),
        width=1600,
        height=400,
    )
    path = diagnose_path(
        False, out, ncount=ncount, neutrons=neutrons, seed=seed, transmission=transmission
    )
    return [top, *path]


if __name__ == "__main__":
    main(out=parse_out(__doc__))
