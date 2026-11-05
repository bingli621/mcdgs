"""How the bender passes each band of wavelengths, from a scan of its rotation.

The scan runs the :class:`.Extraction` and reads its beam monitor, after guide 4: its time and
wavelength map, summed over time, in bands of :data:`BAND_BINS` wavelength bins.
:func:`run_bender_scan` scans ``b_rot`` with the bender in, and runs once with it out;
:func:`plot_bender_scan` draws the intensity per band of the scan over that of the same band
with the bender out, from :func:`transmission`. :func:`read_angles` and :func:`read_spectra`
read the files of the scan. :func:`main` does both; run the module for it::

    pixi run python -m mcdgs.diagnostics.trex.bender_scan
"""

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from mcstastox.LoadFile import Data
from numpy.typing import NDArray

from mcdgs.diagnostics.cli import output_folder, parse_out
from mcdgs.diagnostics.trex.extraction import Extraction
from mcdgs.instruments.trex.guides.bispectral import BENDER_ROTATION
from mcdgs.mcstas.simulation import Simulation

BAND_BINS = 8  # wavelength bins of the beam monitor per band: 12 bands of 0.442 Å, 0 to 5.3 Å


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
    ratio = intensity / reference
    ratio_error = ratio * np.sqrt((error / intensity) ** 2 + (reference_error / reference) ** 2)
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
    :data:`BAND_BINS` wavelength bins. The bins are independent, so their errors add in
    quadrature.

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
        If the wavelength bins of the beam monitor are not a whole number of bands.
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

    if len(wavelengths) % BAND_BINS:
        raise ValueError(
            f"{len(wavelengths)} wavelength bins are not a whole number of bands of {BAND_BINS}"
        )
    bands = (-1, BAND_BINS)  # the last axis, in bands of BAND_BINS bins
    intensity = np.array(intensities).reshape(len(intensities), *bands).sum(axis=-1)
    error = np.sqrt(np.array(variances).reshape(len(variances), *bands).sum(axis=-1))
    return wavelengths.reshape(bands).mean(axis=-1), intensity, error


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

    simulation = Simulation(Extraction(), work_dir=work_dir)
    scan_file = simulation.scan(
        "b_rot", angles, "bender_scan", ncount=ncount, seed=seed, bender_in=1
    )
    print(f"Wrote {scan_file}")
    bender_out_file = simulation.run("bender_out", ncount=ncount, seed=seed, bender_in=0)
    print(f"Wrote {bender_out_file}")
    return scan_file, bender_out_file


def plot_bender_scan(scan_file: Path, bender_out_file: Path) -> Figure:
    """Draw the scan of :func:`run_bender_scan`.

    Each curve is one wavelength band, divided by the same band without the bender, so the
    bender-out intensity is 1.

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
    half_width = (centers[1] - centers[0]) / 2

    fig, ax = plt.subplots()
    for band, center in enumerate(centers):
        label = f"{center:.4g}±{half_width:.4g} Å"
        ax.errorbar(angles, ratio[:, band], ratio_error[:, band], color=f"C{band}", label=label)
    ax.axhline(1, color="black", linestyle="--", label="bender out")
    ax.axvline(
        BENDER_ROTATION, color="black", linestyle=":", label=f"b_rot = {BENDER_ROTATION:.3f} deg"
    )
    ax.set_ylim(bottom=0)
    ax.set_xlim(right=0.8)
    ax.set_xlabel("b_rot [deg]")
    ax.set_ylabel("Intensity at the beam monitor / bender out")
    ax.legend()
    return fig


def main(out: Path | None = None, ncount: float = 1e6, seed: int = 1) -> Path:
    """Run the scan of the bender and save its figure.

    Parameters
    ----------
    out:
        Folder of the figure, ``bender_scan.png``. If ``None``, the folder of the scan, next to
        its ``mccode.h5``.
    ncount:
        Number of McStas neutrons, per angle and in the run without the bender.
    seed:
        Seed of the random numbers.

    Returns
    -------
    :
        The path of the figure.
    """
    scan_file, bender_out_file = run_bender_scan(ncount=ncount, seed=seed)
    figure = output_folder(out, scan_file.parent) / "bender_scan.png"
    fig = plot_bender_scan(scan_file, bender_out_file)
    fig.savefig(figure, dpi=150)
    plt.close(fig)
    print(f"Saved {figure}")
    return figure


if __name__ == "__main__":
    main(out=parse_out(__doc__))
