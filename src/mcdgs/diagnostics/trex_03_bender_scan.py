"""How the bender passes each band of wavelengths, from a scan of its rotation.

The scan runs the :class:`.Extraction` and reads its beam monitor, after guide 4: its time and
wavelength map, summed over time, in bands of :data:`BAND_WIDTH` centered on its multiples.
:func:`run_bender_scan` scans ``b_rot`` with the bender in, and runs once with it out;
:func:`plot_bender_scan` draws the intensity per band of the scan over that of the same band
with the bender out, from :func:`transmission`. :func:`read_angles` and :func:`read_spectra`
read the files of the scan. :func:`main` does both, and runs the path with the bender in at its
nominal rotation with :func:`.diagnose_path`; run the module for it::

    pixi run python -m mcdgs.diagnostics.trex_03_bender_scan
"""

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from mcstastox.LoadFile import Data
from numpy.typing import NDArray

from mcdgs.diagnostics.cli import output_folder, parse_out, save_figure
from mcdgs.diagnostics.trex_02_extraction import Extraction, diagnose_path
from mcdgs.instruments.trex.guides.bispectral import BENDER_ROTATION
from mcdgs.mcstas.simulation import Simulation

BAND_WIDTH = 0.5  # Å: the bands are centered on its multiples, 0.5, 1.0, 1.5 Å, ...


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
    work_dir: str | Path | None = None,
) -> tuple[Path, Path]:
    """Scan the rotation of the bender, and run once without it.

    The scan goes to the run folder ``bender_scan`` of the output folder, the run without the
    bender to ``bender_out``; each replaces the one of an earlier call.

    Parameters
    ----------
    angles:
        The values of ``b_rot``, in deg. If ``None``, 19 from -0.8 to 0.1 deg.
    ncount:
        Number of McStas neutrons, per angle and in the run without the bender.
    seed:
        Seed of the random numbers.
    work_dir:
        Work folder, of the input and the output folders. If ``None``, the current folder.

    Returns
    -------
    :
        The files of the scan and of the run without the bender.
    """
    if angles is None:
        angles = np.linspace(-0.8, 0.1, 19)

    # --8<-- [start:run_scan]
    simulation = Simulation(Extraction(), work_dir=work_dir)
    scan_file = simulation.scan(
        "b_rot", angles, "bender_scan", ncount=ncount, seed=seed, bender_in=1
    )
    print(f"Saved {scan_file}")
    bender_out_file = simulation.run("bender_out", ncount=ncount, seed=seed, bender_in=0)
    # --8<-- [end:run_scan]
    print(f"Saved {bender_out_file}")
    return scan_file, bender_out_file


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


def main(
    out: Path | None = None,
    ncount: float = 1e6,
    seed: int = 1,
    path_ncount: float = 1e7,
    neutrons: int = 1_000_000,
    transmission: float = 0.2,
) -> list[Path]:
    """Run the scan of the bender and save its figure, then run the path with the bender in.

    Parameters
    ----------
    out:
        Folder of the figures. If ``None``, the folder of the scan, next to its ``mccode.h5``,
        for ``bender_scan.png``, and that of the McStas run of the path for the others.
    ncount:
        Number of McStas neutrons, per angle and in the run without the bender.
    seed:
        Seed of the random numbers.
    path_ncount:
        Number of McStas neutrons in each run of the path with the bender in, at its nominal
        rotation: see :func:`.diagnose_path`.
    neutrons, transmission:
        Number of tof neutrons, and the factor on tof at the beam monitor, of that path.

    Returns
    -------
    :
        The path of the figure of the scan, then those of the diagnostics, the diagram and the
        comparison with the bender in.
    """
    scan_file, bender_out_file = run_bender_scan(ncount=ncount, seed=seed)
    fig = plot_bender_scan(scan_file, bender_out_file)
    figure = save_figure(fig, output_folder(out, scan_file.parent) / "bender_scan.png")
    # --8<-- [start:bender_in]
    path = diagnose_path(
        True, out, ncount=path_ncount, neutrons=neutrons, seed=seed, transmission=transmission
    )
    # --8<-- [end:bender_in]
    return [figure, *path]


if __name__ == "__main__":
    main(out=parse_out(__doc__))
