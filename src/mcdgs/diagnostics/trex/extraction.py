"""T-Rex up to guide 4, the thermal and the cold path, to study the extraction guides.

:class:`Extraction` is the section from the source to the beam monitor after guide 4, at 5.9 m,
with the same components as :class:`.TRex`, but for the band: its source and its beam monitor
keep the widest band, 0.1 to 10 Å, where those of :class:`.TRex` follow ``lambda0``. It also has
the source monitor of :class:`.TRexSource`, at the entrance of guide 1.

:func:`main` compares McStas and tof at the beam monitor, with the bender out and in. tof has
no guides: every neutron of its pulse reaches the monitor. tof is scaled to the total of McStas
at the source monitor, where both see the same beam, and by ``transmission`` on top of it, so
the beam monitor shows what the guides let through in McStas against that factor. The McStas runs also draw McStasScript's diagnostics of each path: the beam
after the source, after guide 3 and at the exit of guide 4b. Run the module for both::

    pixi run python -m mcdgs.diagnostics.trex.extraction
"""

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import mcstasscript as ms
from matplotlib.figure import Figure

from mcdgs.core.base import Instrument
from mcdgs.core.components import (
    Arm,
    ESS_butterfly,
    Guide_gravity,
    Monitor_nD,
    Pol_bender,
    Progress_bar,
)
from mcdgs.core.params import InstrParam, param
from mcdgs.diagnostics.cli import output_folder, parse_out
from mcdgs.diagnostics.compare import plot_comparison, read_comparison
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
from mcdgs.mcstas.builder import build_mcstas, default_input_dir
from mcdgs.mcstas.run_folders import default_output_dir, run_folder
from mcdgs.tof.builder import build_tof
from mcdgs.tof.reader import monitor_bins


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
        for name, component in guides_extraction().items():
            setattr(self, name, component)
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


def _save(fig: Figure, path: Path, dpi: int) -> Path:
    """Save a figure, close it and return its path."""
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    print(f"Saved {path}")
    return path


def main(
    out: Path | None = None,
    ncount: float = 1e6,
    neutrons: int = 1_000_000,
    seed: int = 1,
    transmission: float = 1.0,
) -> list[Path]:
    """Compare McStas and tof at the beam monitor, bender out and in, and save the figures.

    McStas runs once per path, with McStasScript's diagnostics, and tof once: tof has no bender,
    so the same run serves both. On each path, tof is scaled to the total of McStas at the source
    monitor, and by ``transmission`` on top of it. Both figures of a
    path, the diagnostics and the comparison,
    ``trex_extraction_bender_out_diagnose.png`` and ``trex_extraction_bender_out.png``, or
    ``..._in``, go to ``out``.

    Parameters
    ----------
    out:
        Folder of the figures. If ``None``, the folder of the McStas run of each path.
    ncount:
        Number of McStas neutrons, in each run.
    neutrons:
        Number of tof neutrons, in each run.
    seed:
        Seed of the random numbers of both models.
    transmission:
        The factor on tof at the beam monitor, on top of the scale of the source monitor: the
        fraction of the beam the guides would let through, which tof does not model.

    Returns
    -------
    :
        The paths of the figures, bender out then in: the diagnostics, then the comparison.
    """
    extraction = Extraction()

    # the bin edges of each monitor, read from it so that they cannot drift apart
    source_bins = monitor_bins(extraction.source_monitor)
    bins = monitor_bins(extraction.beam_monitor)

    tof_result = build_tof(extraction, neutrons=neutrons, seed=seed).run()

    # McStasScript's diagnostics, the same for both paths: the beam after the source, guide 3
    # and guide 4b, its profile, wavelength against x, and its spectrum
    diag = ms.Diagnostics(build_mcstas(extraction, input_path=default_input_dir(extraction)))
    for component in ("source", "guide_3", "guide_4b"):
        # rays="all": a limited list stops early, which leaves the intensities too low
        diag.add_point(after=component, rays="all")  # pyright: ignore[reportArgumentType]
    diag.add_view("x", "y", bins=100)
    diag.add_view("x", "l", bins=100)
    diag.add_view("l")

    figures = []
    for bender_in, path in ((False, "bender out, guide 4a"), (True, "bender in, guide 4b")):
        run_name = f"trex_extraction_bender_{'in' if bender_in else 'out'}"
        folder = run_folder(default_output_dir(), run_name, "diagnose")
        diag.set_parameters(bender_in=int(bender_in))
        diag.settings(
            ncount=ncount, suppress_output=True, mpi="auto", seed=seed, output_path=str(folder)
        )
        diag.run()  # a whole McStas run: its folder also holds the beam monitor
        diag.plot()  # draws on a new figure, which it does not return
        figure_folder = output_folder(out, folder)
        figures.append(_save(plt.gcf(), figure_folder / f"{run_name}_diagnose.png", dpi=100))

        # tof's counts are scaled to the total of McStas at the source monitor, in n/s, and by
        # the transmission on top of it at the beam monitor
        mcstas_data = ms.load_data(str(folder))
        at_source = read_comparison(
            "source_monitor", mcstas_data, tof_result, source_bins, ncount, neutrons
        )
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
        mcstas_total = comparison.mcstas_intensity.sum()
        print(
            f"{path}: McStas total {mcstas_total:.4g} n/s, "
            f"{mcstas_total / at_source.mcstas_intensity.sum():.4g} of that at the source, "
            f"{comparison.scale:.4g} n/s per tof neutron"
        )
        fig = plot_comparison(comparison, f"T-Rex flux at the beam monitor after guide 4, {path}")
        figures.append(_save(fig, figure_folder / f"{run_name}.png", dpi=150))
    return figures


if __name__ == "__main__":
    main(out=parse_out(__doc__))
