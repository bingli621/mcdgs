"""Scan the rotation of the bender and see how it passes each band of wavelengths.

The beam is measured at the exit of guide 4b, with one monitor bin per band. A run without the
bender, through the thermal guide 4a, gives the intensity to compare with. The figure has one
curve per band, divided by that intensity.

:func:`simulate` runs McStas into two folders of mcstas_output/, replacing the earlier ones.
:func:`plot` only reads them, and draws the figure of the documentation. To replot without
simulating, leave out ``simulate()`` at the end of this file.
"""

from dataclasses import dataclass
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from mcstastox.LoadFile import Data
from numpy.typing import NDArray

from mcdgs.core.base import Instrument
from mcdgs.core.components import Guide_gravity, Monitor_nD, Progress_bar
from mcdgs.core.params import InstrParam, param
from mcdgs.instruments.trex.guides.extraction import (
    BENDER_ROTATION,
    bender,
    guides_extraction,
)
from mcdgs.instruments.trex.sources import source
from mcdgs.mcstas.run_folders import default_output_dir
from mcdgs.mcstas.simulation import Simulation

WORK_DIR = Path(__file__).resolve().parents[2]  # the repository, wherever the script runs from
FIGURE = WORK_DIR / "docs" / "instruments" / "trex" / "figures" / "diagnose_bender.png"
SCAN_FILE = default_output_dir(WORK_DIR) / "bender_scan" / "mccode.h5"
BENDER_OUT_FILE = default_output_dir(WORK_DIR) / "bender_out" / "mccode.h5"
WAVELENGTH_EDGES = np.arange(0.5, 5.25, 0.4)  # Å, evenly spaced: the bins of the monitor


@dataclass(frozen=True)
class Params:
    """Control parameters of the scan.

    Frozen: change them with ``BenderScan.update_params``. Each field's unit and comment are in
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


DEFAULT_PARAMS = Params()  # frozen, so one shared default is safe


class BenderScan(Instrument):
    """T-Rex up to guide 4, with a wavelength spectrum at the exit of guide 4b.

    Parameters
    ----------
    params:
        Control parameters.
    """

    params_class = Params
    name = "trex_bender_scan"  # the McStas instrument name

    guide_4b: Guide_gravity  # set by name below, from guides_extraction

    def __init__(self, params: Params = DEFAULT_PARAMS) -> None:
        super().__init__(params)
        self.origin = Progress_bar(at=[0, 0, 0])
        self.source = source()
        to_guide_4a, to_guide_4b = guides_extraction()
        for name, component in {**to_guide_4a, "bender": bender(), **to_guide_4b}.items():
            setattr(self, name, component)
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

    @property
    def instr_params(self) -> list[InstrParam]:
        return [
            InstrParam("bender_in", int(self.params.bender_in), "", "bender in = 1, out = 0"),
            InstrParam("b_rot", self.params.b_rot, "deg", "rotation of the cold neutron bender"),
        ]


def simulate() -> None:
    """Scan the rotation of the bender and run once without it."""
    angles = np.linspace(-0.8, 0.1, 19)  # deg, values of b_rot
    ncount = 1e7  # neutrons per angle
    seed = 1

    simulation = Simulation(BenderScan(), work_dir=WORK_DIR)
    folder_name = SCAN_FILE.parent.name
    simulation.scan("b_rot", angles, folder_name, ncount=ncount, seed=seed, bender_in=1)
    print(f"Wrote {SCAN_FILE}")
    folder_name = BENDER_OUT_FILE.parent.name
    simulation.run(folder_name, ncount=ncount, seed=seed, bender_in=0)
    print(f"Wrote {BENDER_OUT_FILE}")


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
    """Read the intensity per wavelength band of every scan point.

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


def plot() -> None:
    """Plot the intensity against the rotation of the bender, one curve per wavelength band.

    Each curve is divided by the intensity of its band without the bender, which is therefore 1.
    Everything comes from the files that :func:`simulate` writes.
    """
    angles = read_angles(SCAN_FILE)
    centers, intensity, error = read_spectra(SCAN_FILE)
    _, bender_out, bender_out_error = read_spectra(BENDER_OUT_FILE)
    reference, reference_error = bender_out[0], bender_out_error[0]

    ratio = intensity / reference
    ratio_error = ratio * np.sqrt((error / intensity) ** 2 + (reference_error / reference) ** 2)
    half_width = (centers[1] - centers[0]) / 2

    fig, ax = plt.subplots()
    for band, center in enumerate(centers):
        label = f"{center:.4g}±{half_width:.4g} Å"
        ax.errorbar(angles, ratio[:, band], ratio_error[:, band], color=f"C{band}", label=label)
    ax.axhline(1, color="black", linestyle="--", label="bender out")
    ax.axvline(
        BENDER_ROTATION,
        color="black",
        linestyle=":",
        label=f"b_rot = {BENDER_ROTATION:.3f} deg",
    )
    ax.set_ylim(bottom=0)
    ax.set_xlabel("b_rot [deg]")
    ax.set_ylabel("Intensity after guide 4b / bender out")
    ax.legend()
    # ax.spines[["top", "right"]].set_visible(False)
    # fig.legend(loc="outside right upper", frameon=False)
    FIGURE.parent.mkdir(exist_ok=True)
    fig.savefig(FIGURE, dpi=150)
    print(f"Saved {FIGURE}")


if __name__ == "__main__":
    # simulate()
    plot()
    plt.show()
