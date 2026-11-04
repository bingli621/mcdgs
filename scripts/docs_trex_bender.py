"""Make the figures of the T-Rex extraction page, from McStas runs of T-Rex up to guide 4.

`pixi run docs-bender` runs this, and `pixi run docs-bender --simulate` runs the scan again
too. Run it when the extraction guides or the bender change, and commit what it writes.

- docs/instruments/trex/figures/diag_after_4b_bender_out.png and ..._bender_in.png: the beam
  after the source, after guide 3 and at the exit of guide 4, for each path. McStasScript's
  diagnostics run McStas for these every time, with 1e6 neutrons each.
- docs/instruments/trex/figures/diagnose_bender.png: the intensity after guide 4b against the
  rotation of the bender, one curve per wavelength band, divided by the intensity with the
  bender out. It is drawn from the last scan in mcstas_output/: the scan, 19 angles of 1e7
  neutrons and one run with the bender out, is long, so it runs only with --simulate, and
  replaces the earlier one.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # no window: the figures only go to files

import matplotlib.pyplot as plt
import mcstasscript as ms
import numpy as np

from mcdgs.diagnostics.trex.bender import BenderScan, read_angles, read_spectra, transmission
from mcdgs.diagnostics.trex.extraction import Extraction, ExtractionParams
from mcdgs.instruments.trex.guides.bispectral import BENDER_ROTATION
from mcdgs.mcstas.builder import build_mcstas, default_input_dir
from mcdgs.mcstas.run_folders import default_output_dir, run_folder
from mcdgs.mcstas.simulation import Simulation

REPO = Path(__file__).resolve().parents[1]
FIGURES = REPO / "docs" / "instruments" / "trex" / "figures"
SCAN_FIGURE = FIGURES / "diagnose_bender.png"
SCAN_FILE = default_output_dir(REPO) / "bender_scan" / "mccode.h5"
BENDER_OUT_FILE = default_output_dir(REPO) / "bender_out" / "mccode.h5"
SEED = 1


def diagnose(bender_in: bool) -> None:
    """Run the diagnostics of one path and save its figure.

    The beam after the source, after guide 3 and at the exit of guide 4b: its profile,
    wavelength against x, and its spectrum.

    Parameters
    ----------
    bender_in:
        True for the cold path, through the bender and guide 4b; False for the thermal one,
        through guide 4a.
    """
    extraction = Extraction(ExtractionParams(bender_in=bender_in))
    mcstas_instr = build_mcstas(extraction, input_path=default_input_dir(extraction, REPO))
    mcstas_instr.set_parameters(bender_in=int(bender_in))

    diag = ms.Diagnostics(mcstas_instr)
    output_path = run_folder(default_output_dir(REPO), None, "diagnose")
    diag.settings(
        ncount=1e6, suppress_output=True, mpi="auto", seed=SEED, output_path=str(output_path)
    )
    for component in ("source", "guide_3", "guide_4b"):
        # rays="all": a limited list stops early, which leaves the intensities too low
        diag.add_point(after=component, rays="all")  # pyright: ignore[reportArgumentType]
    diag.run()

    diag.add_view("x", "y", bins=100)
    diag.add_view("x", "l", bins=100)
    diag.add_view("l")
    diag.plot()  # draws on a new figure, which it does not return
    figure = FIGURES / f"diag_after_4b_bender_{'in' if bender_in else 'out'}.png"
    plt.gcf().savefig(figure, dpi=100)
    plt.close("all")
    print(f"Saved {figure}")


def simulate() -> None:
    """Scan the rotation of the bender, and run once without it, into mcstas_output/."""
    angles = np.linspace(-0.8, 0.1, 19)  # deg, values of b_rot
    ncount = 1e7  # neutrons per angle

    simulation = Simulation(BenderScan(), work_dir=REPO)
    simulation.scan("b_rot", angles, SCAN_FILE.parent.name, ncount=ncount, seed=SEED, bender_in=1)
    print(f"Wrote {SCAN_FILE}")
    simulation.run(BENDER_OUT_FILE.parent.name, ncount=ncount, seed=SEED, bender_in=0)
    print(f"Wrote {BENDER_OUT_FILE}")


def plot_scan() -> None:
    """Plot the scan of :func:`simulate` and save the figure.

    Each curve is one wavelength band, divided by the same band without the bender, so the
    bender-out intensity is 1.
    """
    angles = read_angles(SCAN_FILE)
    centers, intensity, error = read_spectra(SCAN_FILE)
    _, bender_out, bender_out_error = read_spectra(BENDER_OUT_FILE)
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
    ax.set_xlabel("b_rot [deg]")
    ax.set_ylabel("Intensity after guide 4b / bender out")
    ax.legend()
    fig.savefig(SCAN_FIGURE, dpi=150)
    plt.close(fig)
    print(f"Saved {SCAN_FIGURE}")


def main() -> None:
    """Run the diagnostics of both paths, the scan if asked, and draw every figure."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--simulate", action="store_true", help="run the scan again: long, replaces the last"
    )
    args = parser.parse_args()

    FIGURES.mkdir(parents=True, exist_ok=True)
    for bender_in in (False, True):
        diagnose(bender_in)
    if args.simulate:
        simulate()
    plot_scan()


if __name__ == "__main__":
    main()
