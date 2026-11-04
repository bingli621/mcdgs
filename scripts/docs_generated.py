"""Make the figures and results of the documentation that need only mcdgs and tof.

`pixi run -e docs docs-build` and `docs` run this first, so these are always made from the
current code. They go to docs/instruments/trex/generated/, which git ignores. The figures and
files that need McStas are made by the notebooks instead, and committed.

- trex_print.txt: print(trex), the parameters and the position of every component;
- tof_result.txt: print() of the tof run at the default control parameters;
- tof_diagram.png: its time-distance diagram;
- tof_modes.txt: what reaches the sample in each mode, as a Markdown table.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import get_args

import matplotlib

matplotlib.use("Agg")  # no window: the figure only goes to a file

import numpy as np
import tof

from mcdgs.instruments.trex.choppers import Mode
from mcdgs.instruments.trex.model import TRex
from mcdgs.tof.builder import build_tof

OUT = Path(__file__).resolve().parents[1] / "docs" / "instruments" / "trex" / "generated"

NEUTRONS = 1_000_000
SEED = 1  # the same neutrons in every build, so the results and the pages agree
SAMPLE = "sample_pos"

# sub-pulses are told apart by the gaps between arrivals at the sample: at the default
# parameters these are about 4 ms (one turn of the M-Chopper) between sub-pulses and under
# 0.01 ms within one
SUB_PULSE_GAP_MS = 1.0


@dataclass(frozen=True)
class AtSample:
    """What reaches the sample in one tof run."""

    neutrons: int
    wavelength_min: float  # Å
    wavelength_max: float  # Å
    sub_pulses: int
    first_arrival: float  # ms after the start of the pulse
    last_arrival: float  # ms after the start of the pulse


def run(trex: TRex) -> tof.Result:
    """Run tof on T-Rex, with the sample position as the detector.

    Parameters
    ----------
    trex:
        The instrument, at the control parameters to run.

    Returns
    -------
    :
        The result of the run.
    """
    return build_tof(trex, neutrons=NEUTRONS, detectors=(SAMPLE,), seed=SEED).run()


def at_sample(result: tof.Result) -> AtSample:
    """Summarize the neutrons of a run that reach the sample.

    Parameters
    ----------
    result:
        The result of :func:`run`.

    Returns
    -------
    :
        Their number, wavelength band, sub-pulses and arrival times.
    """
    neutrons = result[SAMPLE].data["pulse", 0]
    reached = neutrons[~neutrons.masks["blocked_by_others"]]
    wavelengths = reached.coords["wavelength"].to(unit="angstrom").values
    arrivals = np.sort(reached.coords["toa"].to(unit="ms").values)
    gaps = np.diff(arrivals)
    return AtSample(
        neutrons=len(arrivals),
        wavelength_min=float(wavelengths.min()),
        wavelength_max=float(wavelengths.max()),
        sub_pulses=int(np.count_nonzero(gaps > SUB_PULSE_GAP_MS)) + 1,
        first_arrival=float(arrivals[0]),
        last_arrival=float(arrivals[-1]),
    )


def modes_table(by_mode: dict[str, AtSample]) -> str:
    """Return a Markdown table of what reaches the sample in each mode.

    Parameters
    ----------
    by_mode:
        The summary of each mode, by mode.

    Returns
    -------
    :
        The table, with a line per mode.
    """
    header = f"Mode | Neutrons at the sample, of {NEUTRONS:,} | Wavelengths, Å | Sub-pulses | Arrival, ms"
    lines = [f"| {header} |", "| ---- | ---- | ---- | ---- | ---- |"]
    for mode, s in by_mode.items():
        lines.append(
            f'| `"{mode}"` | {s.neutrons} | {s.wavelength_min:.2f} to {s.wavelength_max:.2f} '
            f"| {s.sub_pulses} | {s.first_arrival:.0f} to {s.last_arrival:.0f} |"
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    """Write every generated file."""
    OUT.mkdir(parents=True, exist_ok=True)
    trex = TRex()  # the default control parameters
    (OUT / "trex_print.txt").write_text(f"{trex}\n")

    default_mode = trex.params.mode
    result = run(trex)
    (OUT / "tof_result.txt").write_text(f"{str(result).rstrip()}\n")
    diagram = result.plot(visible_rays=1000, blocked_rays=5000)
    diagram.fig.savefig(OUT / "tof_diagram.png", dpi=150)

    by_mode = {default_mode: at_sample(result)}
    for mode in get_args(Mode):
        if mode != default_mode:
            trex.update_params(mode=mode)
            by_mode[mode] = at_sample(run(trex))
    (OUT / "tof_modes.txt").write_text(modes_table(by_mode))
    print(f"Written to {OUT}")


if __name__ == "__main__":
    main()
