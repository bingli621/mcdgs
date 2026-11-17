"""T-Rex up to guide 4, the thermal and the cold path, to study the extraction guides.

:class:`Extraction` is the section from the source to the beam monitor after guide 4, at 5.9 m,
with the same components as :class:`.TRex`, but for the band: its source and its beam monitor
keep the widest band, 0.1 to 10 Å, where those of :class:`.TRex` follow ``lambda0``. It also has
the source monitor of :class:`.TRexSource`, at the entrance of guide 1.

:func:`diagnose_mcstas` runs McStas and :func:`diagnose_tof` compares tof with it at the beam
monitor, on one path through guide 4. tof has no guides: every neutron of its pulse reaches the
monitor. tof is scaled to the total of McStas at the source monitor, where both see the same
beam, and by ``transmission`` on top of it, so the beam monitor shows what the guides let
through in McStas against that factor. The McStas runs also draw McStasScript's diagnostics, the
beam after the source, after guide 3 and at the exit of guide 4b, and its diagram of the
intensity before each component. :func:`main` draws the extraction from above and runs both
paths: without the bender, through guide 4a, and with it, through the bender and guide 4b. Run the module for it::

    pixi run python -m mcdgs.diagnostics.trex.extraction
"""

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import mcstasscript as ms

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
from mcdgs.diagnostics.compare import (
    MonitorComparison,
    compare_maps,
    plot_comparison,
    read_maps,
    tof_scale,
)
from mcdgs.diagnostics.drawing import draw_instrument, draw_intensity_diagram
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
from mcdgs.mcstas import Simulation, build_mcstas
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


# --8<-- [start:instrument]
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


def define_instrument(bender_in: bool = False) -> Extraction:
    """Return the instrument of this diagnostic.

    Parameters
    ----------
    bender_in:
        The path through guide 4: True for the bender and guide 4b, False for guide 4a.

    Returns
    -------
    :
        :class:`Extraction` on that path, at the nominal rotation of the bender.
    """
    params = ExtractionParams(bender_in=bender_in)
    return Extraction(params)


# --8<-- [end:instrument]


def diagnose_mcstas(
    instrument: Extraction,
    run_name: str,
    figure_folder: Path,
    ncount: float = 1e7,
    seed: int = 1,
    diagnostics_ncount: float = DIAGNOSTICS_NCOUNT,
) -> tuple[list[Path], Path]:
    """Run McStas twice on one path through guide 4, and save its two figures.

    The first run has McStasScript's diagnostics, the beam after the source, after guide 3 and
    at the exit of guide 4b; the second its diagram of the intensity before each component,
    which also gives the monitors. The diagnostics list every ray at their three positions, so
    their run has fewer neutrons. The figures are ``<run_name>_diagnose.png`` and
    ``<run_name>_intensity.png``; the runs go to the run folders ``<run_name>`` and
    ``<run_name>_intensity``.

    Parameters
    ----------
    instrument:
        :class:`Extraction` on the path, which the runs take as their default.
    run_name:
        Name of the runs and of the figures.
    figure_folder:
        Folder of the figures.
    ncount:
        Number of McStas neutrons of the run of the diagram and the monitors.
    seed:
        Seed of the random numbers.
    diagnostics_ncount:
        Number of McStas neutrons of the run with the diagnostics.

    Returns
    -------
    :
        The paths of the diagnostics and the diagram, and the folder of the run of the
        monitors.
    """
    # --8<-- [start:diagnostics]
    # McStasScript's diagnostics: the beam after the source, guide 3 and guide 4b, its profile,
    # wavelength against x, and its spectrum. McStasScript runs them itself, so the simulation
    # only prepares the run: its folder, its settings and its record, in McStas's own format,
    # which McStasScript reads
    simulation = Simulation(instrument)
    run, mcstas_instr = simulation.prepare(
        run_name, ncount=diagnostics_ncount, seed=seed, nexus=False, mpi="auto"
    )
    diag = ms.Diagnostics(mcstas_instr)  # runs a copy of mcstas_instr, with its settings
    for component in ("source", "guide_3", "guide_4b"):
        # rays="all": a limited list stops early, which leaves the intensities too low
        diag.add_point(after=component, rays="all")  # pyright: ignore[reportArgumentType]
    diag.add_view("x", "y", bins=100, colormap="turbo")  # McStasScript's default is jet
    diag.add_view("x", "l", bins=100, colormap="turbo")
    diag.add_view("l")
    # --8<-- [end:diagnostics]

    diag.run()
    run.finish()  # the record of the run, in its folder
    diag.plot()  # draws on a new figure, which it does not return
    figures = [save_figure(plt.gcf(), figure_folder / f"{run_name}_diagnose.png", dpi=100)]

    # McStasScript's diagram of the components, with the intensity and the rays before each:
    # a McStas run of its own, in a folder of its own, which also holds the monitors
    run, mcstas_instr = simulation.prepare(
        f"{run_name}_intensity", ncount=ncount, seed=seed, nexus=False, mpi="auto"
    )
    # --8<-- [start:diagram]
    diagram = draw_intensity_diagram(
        mcstas_instr, figure_folder / f"{run_name}_intensity.png", DIAGRAM_INTENSITY_LIMITS
    )
    # --8<-- [end:diagram]
    run = run.finish()
    figures.append(diagram)
    return figures, run.folder


def diagnose_tof(
    instrument: Extraction,
    mcstas_folder: Path,
    ncount: float = 1e7,
    neutrons: int = 1_000_000,
    seed: int = 1,
    transmission: float = 0.2,
) -> tuple[MonitorComparison, MonitorComparison]:
    """Run tof on one path through guide 4, and scale it to McStas.

    tof has no guides and no bender: every neutron of its pulse reaches the beam monitor, and
    its run is the same on both paths. It is scaled to the total of McStas at the source
    monitor, where both see the same beam, and by ``transmission`` on top of it at the beam
    monitor.

    Parameters
    ----------
    instrument:
        :class:`Extraction` on the path.
    mcstas_folder:
        Folder of the McStas run of the monitors, from :func:`diagnose_mcstas`.
    ncount:
        Number of McStas neutrons of the run in ``mcstas_folder``.
    neutrons:
        Number of tof neutrons.
    seed:
        Seed of the random numbers.
    transmission:
        The factor on tof at the beam monitor, on top of the scale of the source monitor: the
        fraction of the beam the guides would let through, which tof does not model. 0.2 is
        an estimate, close to the 0.17 McStas has with the bender out at the defaults, and the
        0.12 with it in.

    Returns
    -------
    :
        Both models at the source monitor, then at the beam monitor.
    """
    tof_result = build_tof(instrument, neutrons=neutrons, seed=seed).run()

    # --8<-- [start:scale_tof]
    # at the source monitor, where both models see the same beam: the n/s each tof neutron
    # stands for, the total of McStas over that of tof
    mcstas_data = ms.load_data(str(mcstas_folder))
    source_bins = monitor_bins(instrument.source_monitor.bin_limits)
    mcstas_source, tof_source = read_maps("source_monitor", mcstas_data, tof_result, source_bins)
    scale = tof_scale(mcstas_source, tof_source)
    at_source = compare_maps("source_monitor", mcstas_source, tof_source, ncount, neutrons, scale)

    # at the beam monitor: tof with that scale, and the transmission of the guides on top of it
    beam_bins = monitor_bins(instrument.beam_monitor.bin_limits)
    mcstas_beam, tof_beam = read_maps("beam_monitor", mcstas_data, tof_result, beam_bins)
    at_beam = compare_maps(
        "beam_monitor", mcstas_beam, tof_beam, ncount, neutrons, scale, transmission
    )
    # --8<-- [end:scale_tof]
    return at_source, at_beam


def main(
    out: Path | None = None,
    ncount: float = 1e7,
    neutrons: int = 1_000_000,
    seed: int = 1,
    transmission: float = 0.2,
) -> list[Path]:
    """Draw the extraction from above, and run and compare both paths through guide 4.

    The view from above is ``trex_extraction_top.png``. Each path, first with the bender out,
    through guide 4a, then with it in, through the bender and guide 4b, runs McStas with
    :func:`diagnose_mcstas` and tof, scaled to it, with :func:`diagnose_tof`, and draws both at
    the beam monitor. The figures of a path are ``trex_extraction_bender_out_diagnose.png``,
    ``trex_extraction_bender_out_intensity.png`` and ``trex_extraction_bender_out.png``, the
    comparison, or ``..._in``.

    Parameters
    ----------
    out:
        Folder of the figures. If ``None``, ``trex_extraction`` in the results folder of
        :meth:`.Workspace.find`.
    ncount, seed:
        As in :func:`diagnose_mcstas`, for both paths; ``seed`` also of tof.
    neutrons, transmission:
        As in :func:`diagnose_tof`, for both paths.

    Returns
    -------
    :
        The paths of the view from above, then of the diagnostics, the diagram and the
        comparison of the path with the bender out, then of those with it in.
    """
    figure_folder = output_folder(out, Extraction.name)
    # from above: z along the beam, x across it, 16 x 4 in, within TOP_VIEW_LIMITS; mcdisplay
    # draws every component, whatever its WHEN
    top = draw_instrument(
        build_mcstas(define_instrument()),
        figure_folder / "trex_extraction_top.png",
        equal_aspect=False,
        projection="zx",
        limits=list(TOP_VIEW_LIMITS),
        width=1600,
        height=400,
    )
    figures = [top]
    # --8<-- [start:paths]
    # each path through guide 4, with the same McStas and tof runs on its own instrument: first
    # through guide 4a, then through the bender and guide 4b
    for bender_in in (False, True):
        instrument = define_instrument(bender_in)
        path = "bender in, guide 4b" if bender_in else "bender out, guide 4a"
        run_name = f"{instrument.name}_bender_{'in' if bender_in else 'out'}"
        mcstas_figures, mcstas_folder = diagnose_mcstas(
            instrument, run_name, figure_folder, ncount=ncount, seed=seed
        )
        at_source, comparison = diagnose_tof(
            instrument,
            mcstas_folder,
            ncount=ncount,
            neutrons=neutrons,
            seed=seed,
            transmission=transmission,
        )

        mcstas_total = comparison.mcstas_intensity.sum()
        print(
            f"{path}: McStas total {mcstas_total:.4g} n/s, "
            f"{mcstas_total / at_source.mcstas_intensity.sum():.4g} of that at the source, "
            f"{comparison.scale:.4g} n/s per tof neutron"
        )
        fig = plot_comparison(comparison, f"T-Rex flux at the beam monitor after guide 4, {path}")
        figures += [*mcstas_figures, save_figure(fig, figure_folder / f"{run_name}.png")]
    # --8<-- [end:paths]
    return figures


if __name__ == "__main__":
    main(out=parse_out(__doc__))
