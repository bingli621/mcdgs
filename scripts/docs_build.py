"""Make the files of the documentation that need only mcdgs and tof.

`pixi run -e docs docs-build` and `docs` run this first (as `docs-generate`), so these are always
made from the current code. They go to docs/instruments/trex/generated/, which git ignores.

- trex_print.txt: print(trex), the parameters and the position of every component;
- tof_result.txt: print() of the tof run at the default control parameters;
- tof_diagram.png: its time-distance diagram;
- tof_modes.txt: what reaches the sample in each mode, as a Markdown table.
"""

from typing import get_args

from docs_common import SEED, TREX_DOCS

from mcdgs.diagnostics.trex.choppers import ArrivalSummary, summarize
from mcdgs.instruments.trex.choppers import Mode
from mcdgs.instruments.trex.model import TRex
from mcdgs.tof.builder import build_tof

OUT = TREX_DOCS / "generated"
NEUTRONS = 1_000_000
SAMPLE = "sample_monitor"


def modes_table(by_mode: dict[str, ArrivalSummary]) -> str:
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
    header = (
        f"Mode | Neutrons at the sample, of {NEUTRONS:,} | Wavelengths, Å | Sub-pulses "
        "| Arrival, ms"
    )
    lines = [f"| {header} |", "| ---- | ---- | ---- | ---- | ---- |"]
    for mode, s in by_mode.items():
        lines.append(
            f'| `"{mode}"` | {s.neutrons} | {s.wavelength_min:.2f} to {s.wavelength_max:.2f} '
            f"| {s.sub_pulses} | {s.first_arrival:.0f} to {s.last_arrival:.0f} |"
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    """Write every file, at the default control parameters."""
    OUT.mkdir(parents=True, exist_ok=True)
    trex = TRex()
    (OUT / "trex_print.txt").write_text(f"{trex}\n")

    default_mode = trex.params.mode
    result = build_tof(trex, neutrons=NEUTRONS, seed=SEED).run()
    (OUT / "tof_result.txt").write_text(f"{str(result).rstrip()}\n")
    diagram = result.plot(visible_rays=1000, blocked_rays=5000)
    diagram.fig.savefig(OUT / "tof_diagram.png", dpi=150)

    # the other modes from the same neutrons, after the default one
    by_mode = {default_mode: summarize(result, SAMPLE)}
    for mode in get_args(Mode):
        if mode != default_mode:
            trex.update_params(mode=mode)
            other = build_tof(trex, neutrons=NEUTRONS, seed=SEED).run()
            by_mode[mode] = summarize(other, SAMPLE)
    (OUT / "tof_modes.txt").write_text(modes_table(by_mode))
    print(f"Written to {OUT}")


if __name__ == "__main__":
    main()
