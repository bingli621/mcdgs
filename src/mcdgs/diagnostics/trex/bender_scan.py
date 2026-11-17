"""How the bender passes each band of wavelengths, from a scan of its rotation.

:class:`TRexBenderScan` is T-Rex up to the beam monitor after guide 4, with the bender in, as the
:class:`.Extraction` has it. The scan reads its beam monitor: its time and wavelength map, summed
over time, in bands of :data:`BAND_WIDTH` centered on its multiples.
:func:`run_bender_scan` scans ``b_rot`` with the bender in, and runs once with it out;
:func:`plot_bender_scan` draws the intensity per band of the scan over that of the same band
with the bender out, from :func:`transmission`. :func:`read_angles` and :func:`read_spectra`
read the files of the scan. :func:`main` does both; run the module for it::

    pixi run python -m mcdgs.diagnostics.trex.bender_scan
"""

from dataclasses import dataclass
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from mcstastox.LoadFile import Data
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
from mcdgs.diagnostics.cli import output_folder, parse_out, save_figure
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
from mcdgs.mcstas import Simulation, Workspace

BAND_WIDTH = 0.5  # Å: the bands are centered on its multiples, 0.5, 1.0, 1.5 Å, ...


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

    The source and the beam monitor keep the widest band, 0.1 to 10 Å. Guide 4a, the bender and
    guide 4b are all in the instrument; ``bender_in`` chooses the path, as an instrument
    parameter, so that a run of the scan can set it.

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


def define_instrument() -> TRexBenderScan:
    """Return the instrument of this diagnostic: :class:`TRexBenderScan`, the bender in.

    The runs of the scan still set ``bender_in`` themselves: the scan with the bender in, the run
    without it with it out.
    """
    params = TRexBenderScanParams()  # the default control parameters: the bender in
    return TRexBenderScan(params)


# --8<-- [end:instrument]


def transmission(
    intensity: NDArray, error: NDArray, reference: NDArray, reference_error: NDArray
) -> tuple[NDArray, NDArray]:
    """Divide the intensity of each band by that of the same band in a reference run.

    The errors are independent, so their relative errors add in quadrature.

    Parameters
    ----------
    intensity:
        The intensity per band of every scan point, with the shape (scan points, bands).
    error:
        Its error, with the same shape.
    reference:
        The intensity per band of the reference run, such as the bender out, with the shape
        (bands,).
    reference_error:
        Its error, with the same shape.

    Returns
    -------
    :
        The ratio and its error, both with the shape of ``intensity``.
    """
    # --8<-- [start:transmission]
    ratio = intensity / reference
    ratio_error = ratio * np.sqrt((error / intensity) ** 2 + (reference_error / reference) ** 2)
    # --8<-- [end:transmission]
    return ratio, ratio_error


def read_angles(h5_file: Path) -> NDArray:
    """Read the rotation of the bender of every scan point.

    Parameters
    ----------
    h5_file:
        The file of the scan.

    Returns
    -------
    :
        The rotations in deg.
    """
    folder = str(h5_file.parent)
    with Data(folder) as data:
        n_entries = data.get_number_of_entries()
    angles = []
    for entry in range(1, n_entries + 1):
        with Data(folder, entry_number=entry) as data:
            angles.append(data.get_instrument_parameters()["b_rot"])
    return np.array(angles, dtype=float)


def read_spectra(h5_file: Path) -> tuple[NDArray, NDArray, NDArray]:
    """Read the intensity per wavelength band at the beam monitor of every scan point.

    The map of the beam monitor is summed over time, then over each band of
    :data:`BAND_WIDTH`, centered on a multiple of it. A wavelength bin goes to the band whose
    center is nearest to its own, and to the band above when it is centered on the edge of two,
    so a band's edges are those of its bins, within half a bin of the band's. Only the bands
    the axis of the monitor covers whole are kept. The bins are independent, so their errors add
    in quadrature.

    Parameters
    ----------
    h5_file:
        The file of the scan.

    Returns
    -------
    :
        The wavelength at the center of each band in Å, then the intensity and its error, both
        with the shape (scan points, bands).

    Raises
    ------
    ValueError
        If the axis of the beam monitor covers no band whole.
    """
    wavelengths, intensities, variances = np.empty(0), [], []
    with h5py.File(h5_file) as file:
        for n in range(1, len(file) + 1):  # by number: by name, entry10 comes before entry2
            components = file[f"entry{n}/instrument/components"]
            assert isinstance(components, h5py.Group)
            monitor = next(c for name, c in components.items() if name.endswith("_beam_monitor"))
            time_map = monitor["output/beam_monitor_t_L"]  # (time, wavelength)
            intensities.append(time_map["data"][()].sum(axis=0))
            variances.append((time_map["errors"][()] ** 2).sum(axis=0))
            wavelengths = time_map["Wavelength__Angs_"][()]  # the bin centers

    half_bin = (wavelengths[1] - wavelengths[0]) / 2
    start, stop = wavelengths[0] - half_bin, wavelengths[-1] + half_bin  # the axis, in Å
    # --8<-- [start:bands]
    # the nearest center, in widths; a bin centered on the edge of two bands, as 0.25 Å, goes to
    # the band above, whatever the rounding of its center
    band_of_bin = np.floor(wavelengths / BAND_WIDTH + 0.5 + 1e-9).astype(int)
    # --8<-- [end:bands]
    tolerance = 1e-9  # Å: an edge of the axis on that of a band covers it
    whole = [
        int(band)
        for band in np.unique(band_of_bin)
        if (band - 0.5) * BAND_WIDTH >= start - tolerance
        and (band + 0.5) * BAND_WIDTH <= stop + tolerance
    ]
    if not whole:
        raise ValueError(
            f"the wavelength axis, {start:.4g} to {stop:.4g} Å, covers no band of {BAND_WIDTH} Å"
        )
    # --8<-- [start:sum_bands]
    intensity_bins, variance_bins = np.array(intensities), np.array(variances)
    intensity = np.stack([intensity_bins[:, band_of_bin == b].sum(axis=1) for b in whole], 1)
    variance = np.stack([variance_bins[:, band_of_bin == b].sum(axis=1) for b in whole], 1)
    # --8<-- [end:sum_bands]
    return np.array(whole) * BAND_WIDTH, intensity, np.sqrt(variance)


def run_bender_scan(
    angles: NDArray | None = None,
    ncount: float = 1e6,
    seed: int = 1,
    workspace: Workspace | None = None,
) -> tuple[Path, Path]:
    """Scan the rotation of the bender, and run once without it.

    The scan goes to the run folder ``bender_scan`` of the workspace, the run without the
    bender to ``bender_out``; each replaces the one of an earlier call.

    Parameters
    ----------
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
        The files of the scan and of the run without the bender.
    """
    if angles is None:
        angles = np.linspace(-0.8, 0.1, 19)

    # --8<-- [start:run_scan]
    sim = Simulation(define_instrument(), workspace)
    scan = sim.scan("b_rot", angles, "bender_scan", ncount=ncount, seed=seed, bender_in=1)
    print(f"Saved {scan.data_file}")
    bender_out = sim.run(folder_name="bender_out", ncount=ncount, seed=seed, bender_in=0)
    # --8<-- [end:run_scan]
    print(f"Saved {bender_out.data_file}")
    return scan.data_file, bender_out.data_file


def plot_bender_scan(scan_file: Path, bender_out_file: Path) -> Figure:
    """Draw the scan of :func:`run_bender_scan`.

    Each curve is one wavelength band, divided by the same band without the bender, so the
    bender-out intensity is 1. The bands take matplotlib's default colors in turn, solid for
    the first ten and dash-dotted for the next, and the legend gives the center and half width
    of each.

    Parameters
    ----------
    scan_file:
        The file of the scan.
    bender_out_file:
        The file of the run without the bender.

    Returns
    -------
    :
        The figure.
    """
    angles = read_angles(scan_file)
    centers, intensity, error = read_spectra(scan_file)
    _, bender_out, bender_out_error = read_spectra(bender_out_file)
    reference, reference_error = bender_out[0], bender_out_error[0]

    ratio, ratio_error = transmission(intensity, error, reference, reference_error)

    fig, ax = plt.subplots(figsize=(8, 4.8), layout="constrained")
    n_colors = 10  # the default color cycle, C0 to C9: past it the colors come again
    for band, center in enumerate(centers):
        ax.errorbar(
            angles,
            ratio[:, band],
            ratio_error[:, band],
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
    scan_file, bender_out_file = run_bender_scan(ncount=ncount, seed=seed)
    fig = plot_bender_scan(scan_file, bender_out_file)
    return save_figure(fig, output_folder(out, "trex_bender_scan") / "bender_scan.png")


if __name__ == "__main__":
    main(out=parse_out(__doc__))
