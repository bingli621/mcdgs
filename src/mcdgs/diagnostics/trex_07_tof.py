"""What the chopper cascade of T-Rex passes to a detector, from a tof run.

The pulse-shaping and monochromatic choppers cut the band of the bandwidth choppers, so the
neutrons arrive in
sub-pulses. :func:`count_groups` tells them apart by the gaps between arrival times, and
:func:`summarize` gives how many neutrons reach a detector, which wavelengths and when, as an
:class:`ArrivalSummary`. :func:`modes_table` puts the summary of each mode in a Markdown table.

:func:`main` writes the files of the documentation that need only mcdgs and tof: the tof run at
the default control parameters, and the table of the modes. Run the module for them::

    pixi run python -m mcdgs.diagnostics.trex_07_tof --out docs/instruments/trex/figures
"""

from dataclasses import dataclass
from pathlib import Path
from typing import get_args

import numpy as np
import tof
from numpy.typing import ArrayLike

from mcdgs.diagnostics.cli import output_folder, parse_out
from mcdgs.instruments.trex.choppers import Mode
from mcdgs.instruments.trex.model import TRex
from mcdgs.tof.builder import build_tof


@dataclass(frozen=True)
class ArrivalSummary:
    """What reaches a detector in one pulse of a tof run.

    Parameters
    ----------
    neutrons:
        How many neutrons reach it.
    wavelength_min, wavelength_max:
        Their shortest and longest wavelength, in Å.
    sub_pulses:
        How many groups they arrive in.
    first_arrival, last_arrival:
        When the first and the last arrive, in ms after the start of the pulse.
    """

    neutrons: int
    wavelength_min: float
    wavelength_max: float
    sub_pulses: int
    first_arrival: float
    last_arrival: float


def count_groups(times: ArrayLike, gap: float) -> int:
    """Count the groups of times, a new one starting after each gap longer than ``gap``.

    Parameters
    ----------
    times:
        The times, in any order. There must be at least one.
    gap:
        The longest gap within a group, in the unit of ``times``.

    Returns
    -------
    :
        The number of groups.

    Examples
    --------
    >>> count_groups([0.0, 0.1, 5.0, 5.2, 10.0], gap=1.0)
    3
    """
    ordered = np.sort(np.asarray(times, dtype=float))
    if ordered.size == 0:
        raise ValueError("there are no times to group")
    return int(np.count_nonzero(np.diff(ordered) > gap)) + 1


def summarize(result: tof.Result, detector: str, gap: float = 1.0) -> ArrivalSummary:
    """Summarize the neutrons of the first pulse that reach a detector.

    Parameters
    ----------
    result:
        The result of a tof run.
    detector:
        The name of the detector.
    gap:
        The longest gap between arrivals within a sub-pulse, in ms. The default is chosen for
        T-Rex: its sub-pulses are about 4 ms apart, one turn of the monochromatic choppers,
        while within one the gaps are under 0.01 ms.

    Returns
    -------
    :
        Their number, wavelength band, sub-pulses and arrival times.
    """
    neutrons = result[detector].data["pulse", 0]
    reached = neutrons[~neutrons.masks["blocked_by_others"]]
    if reached.size == 0:
        raise ValueError(f"no neutron reaches {detector!r}")
    wavelengths = reached.coords["wavelength"].to(unit="angstrom").values
    arrivals = reached.coords["toa"].to(unit="ms").values
    return ArrivalSummary(
        neutrons=int(reached.size),
        wavelength_min=float(wavelengths.min()),
        wavelength_max=float(wavelengths.max()),
        sub_pulses=count_groups(arrivals, gap),
        first_arrival=float(arrivals.min()),
        last_arrival=float(arrivals.max()),
    )


def modes_table(by_mode: dict[str, ArrivalSummary], neutrons: int) -> str:
    """Return a Markdown table of what reaches the sample in each mode.

    Parameters
    ----------
    by_mode:
        The summary of each mode, by mode.
    neutrons:
        Number of tof neutrons of each run, for the header.

    Returns
    -------
    :
        The table, with a line per mode.
    """
    header = (
        f"Mode | Neutrons at the sample, of {neutrons:,} | Wavelengths, Å | Sub-pulses "
        "| Arrival, ms"
    )
    lines = [f"| {header} |", "| ---- | ---- | ---- | ---- | ---- |"]
    for mode, s in by_mode.items():
        lines.append(
            f'| `"{mode}"` | {s.neutrons} | {s.wavelength_min:.2f} to {s.wavelength_max:.2f} '
            f"| {s.sub_pulses} | {s.first_arrival:.0f} to {s.last_arrival:.0f} |"
        )
    return "\n".join(lines) + "\n"


def main(out: Path | None = None, neutrons: int = 1_000_000, seed: int = 1) -> list[Path]:
    """Write the tof run and the table of the modes, at the default parameters.

    - ``tof_result.txt``: ``print()`` of the tof run at the default control parameters;
    - ``tof_modes.txt``: what reaches the sample in each mode, from :func:`modes_table`.

    Parameters
    ----------
    out:
        Folder to write to. If ``None``, the current folder.
    neutrons:
        Number of tof neutrons, in each run.
    seed:
        Seed of the random numbers.

    Returns
    -------
    :
        The paths of the files.
    """
    folder = output_folder(out, Path.cwd())
    trex = TRex()
    default_mode = trex.params.mode
    result = build_tof(trex, neutrons=neutrons, seed=seed).run()
    tof_result = folder / "tof_result.txt"
    tof_result.write_text(f"{str(result).rstrip()}\n")

    # the other modes from the same neutrons, after the default one
    by_mode = {default_mode: summarize(result, "sample_monitor")}
    for mode in get_args(Mode):
        if mode != default_mode:
            trex.update_params(mode=mode)
            other = build_tof(trex, neutrons=neutrons, seed=seed).run()
            by_mode[mode] = summarize(other, "sample_monitor")
    tof_modes = folder / "tof_modes.txt"
    tof_modes.write_text(modes_table(by_mode, neutrons))
    print(f"Written to {folder}")
    return [tof_result, tof_modes]


if __name__ == "__main__":
    main(out=parse_out(__doc__))
