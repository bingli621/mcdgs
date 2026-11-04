"""How the bender passes each band of wavelengths, from a scan of its rotation.

:class:`BenderScan` is the :class:`.Extraction` with a wavelength monitor at the exit of
guide 4b, one bin per band of :data:`WAVELENGTH_EDGES`. A scan of ``b_rot`` with the bender in,
and one run with it out, give the intensity per band; :func:`transmission` divides the first by
the second. :func:`read_angles` and :func:`read_spectra` read the files of the scan.
"""

from pathlib import Path

import h5py
import numpy as np
from mcstastox.LoadFile import Data
from numpy.typing import NDArray

from mcdgs.core.components import Monitor_nD

from .extraction import DEFAULT_PARAMS, Extraction, ExtractionParams

WAVELENGTH_EDGES = np.arange(0.5, 5.25, 0.4)  # Å, evenly spaced: the bins of the monitor


class BenderScan(Extraction):
    """T-Rex up to guide 4, with a wavelength spectrum at the exit of guide 4b.

    Parameters
    ----------
    params:
        Control parameters.
    """

    name = "trex_bender_scan"  # the McStas instrument name

    monitor_after_4b: Monitor_nD

    def __init__(self, params: ExtractionParams = DEFAULT_PARAMS) -> None:
        super().__init__(params)
        self.monitor_after_4b = Monitor_nD(
            xwidth=0.1,  # wider than the guide: every neutron leaving it
            yheight=0.1,
            restore_neutron=1,
            options=(
                f"lambda limits=[{WAVELENGTH_EDGES[0]} {WAVELENGTH_EDGES[-1]}] "
                f"bins={len(WAVELENGTH_EDGES) - 1}"
            ),
            filename="monitor_after_guide_4b",
            at=[0, 0, self.guide_4b.l + 1e-3],  # as the Diag monitors
            relative="guide_4b",
        )


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
    """Read the intensity per wavelength band of every scan point of :class:`BenderScan`.

    Parameters
    ----------
    h5_file:
        The file of the scan.

    Returns
    -------
    :
        The wavelength at the center of each band in Å, then the intensity and its error, both
        with the shape (scan points, bands).
    """
    wavelengths, intensities, errors = np.empty(0), [], []
    with h5py.File(h5_file) as file:
        for n in range(1, len(file) + 1):  # by number: by name, entry10 comes before entry2
            components = file[f"entry{n}/instrument/components"]
            assert isinstance(components, h5py.Group)
            monitor = next(
                c for name, c in components.items() if name.endswith("_monitor_after_4b")
            )
            spectrum = monitor["output/monitor_after_guide_4b_L"]
            intensities.append(spectrum["data"][()])
            errors.append(spectrum["errors"][()])
            wavelengths = spectrum["Wavelength__Angs_"][()]
    return wavelengths, np.array(intensities), np.array(errors)
