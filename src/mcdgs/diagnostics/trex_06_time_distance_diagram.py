"""The time-distance diagram of T-Rex, from its tof model.

:func:`main` writes, at the default control parameters:

- ``trex_time_distance.png``: the time-distance diagram of two pulses, from
  :func:`.time_distance_diagram`;
- ``trex_time_distance_params.txt``: the parameters of that run, from :func:`parameters_table`.

It needs only mcdgs and tof. Run the module for them::

    pixi run python -m mcdgs.diagnostics.trex_06_time_distance_diagram --out FOLDER
"""

from pathlib import Path

import matplotlib.pyplot as plt

from mcdgs.core.components import NXdisk_chopper
from mcdgs.core.utils import evaluate
from mcdgs.diagnostics.cli import output_folder, parse_out
from mcdgs.instruments.trex.model import TRex
from mcdgs.tof.builder import beamline
from mcdgs.tof.time_distance import time_distance_diagram

N_PULSES = 2
"""Pulses of the diagram."""
VISIBLE_RAYS = 2000
"""Most rays drawn per pulse, of the neutrons that reach the sample."""
BLOCKED_RAYS = 5000
"""Most rays drawn per pulse, in gray, of the neutrons a chopper stops."""


def parameters_table(
    trex: TRex,
    *,
    neutrons: int,
    seed: int | None,
    n_pulses: int,
    visible_rays: int,
    blocked_rays: int,
) -> str:
    """Return the parameters of a tof run of T-Rex as Markdown tables.

    Three tables: the tof run, with the band of the source; the control parameters that set the
    choppers; and each chopper, with its distance from the source, frequency and delay.

    Parameters
    ----------
    trex:
        The instrument, at the parameters of the run.
    neutrons:
        Number of tof neutrons per pulse.
    seed:
        Seed of the random numbers.
    n_pulses:
        Number of pulses.
    visible_rays, blocked_rays:
        Most rays drawn per pulse, of the neutrons that reach the sample and that are stopped.

    Returns
    -------
    :
        The tables, separated by blank lines.
    """
    source, after_source = beamline(trex)
    params = trex.numeric_instr_params()
    lmin = float(evaluate(source.Lmin, params))
    lmax = float(evaluate(source.Lmax, params))
    run = [
        "| tof run | Value |",
        "| ------- | ----- |",
        f"| Source | ESS pulse of beam port {source.sector}{source.beamline} at {source.frequency:g} Hz |",
        f"| Band of the source | {lmin:.3f} to {lmax:.3f} Å |",
        f"| `neutrons`, per pulse | {neutrons:,} |",
        f"| `n_pulses` | {n_pulses} |",
        f"| `seed` | {seed} |",
        f"| `visible_rays`, per pulse | {visible_rays:,} |",
        f"| `blocked_rays`, per pulse | {blocked_rays:,} |",
    ]
    p = trex.params
    control = [
        "| Control parameter | Value |",
        "| ----------------- | ----- |",
        f'| `mode` | `"{p.mode}"` |',
        f"| `lambda0` | {p.lambda0:g} Å |",
        f"| `rrm` | {p.rrm} |",
        f"| `t0_offset` | {p.t0_offset * 1e3:g} ms |",
    ]
    choppers = [
        "| Chopper | Distance from the source | Frequency | Delay |",
        "| ------- | ------------------------ | --------- | ----- |",
    ]
    for name, component, distance in after_source:
        if isinstance(component, NXdisk_chopper):
            nu = float(evaluate(component.nu, params))
            delay = float(evaluate(component.delay, params))
            frequency = f"{nu:g}".replace("-", "\u2212")  # a minus sign, as in the docs
            choppers.append(
                f"| `{name}` | {distance:.3f} m | {frequency} Hz | {delay * 1e3:.3f} ms |"
            )
    return "\n\n".join("\n".join(table) for table in (run, control, choppers)) + "\n"


def main(out: Path | None = None, neutrons: int = 1_000_000, seed: int = 1) -> list[Path]:
    """Draw the time-distance diagram of T-Rex and list its parameters, at the defaults.

    Parameters
    ----------
    out:
        Folder to write to. If ``None``, the current folder.
    neutrons:
        Number of tof neutrons per pulse.
    seed:
        Seed of the random numbers.

    Returns
    -------
    :
        The paths of the diagram and of the table of its parameters.
    """
    folder = output_folder(out, Path.cwd())
    path = folder / "trex_time_distance.png"

    # --8<-- [start:time_distance]
    trex = TRex()  # the default control parameters
    figure = time_distance_diagram(
        trex,
        neutrons,
        seed=seed,
        n_pulses=N_PULSES,
        visible_rays=VISIBLE_RAYS,
        blocked_rays=BLOCKED_RAYS,
    )
    figure.savefig(path, dpi=150)
    # --8<-- [end:time_distance]
    plt.close(figure)

    table = folder / "trex_time_distance_params.txt"
    table.write_text(
        parameters_table(
            trex,
            neutrons=neutrons,
            seed=seed,
            n_pulses=N_PULSES,
            visible_rays=VISIBLE_RAYS,
            blocked_rays=BLOCKED_RAYS,
        )
    )
    print(f"Written to {folder}")
    return [path, table]


if __name__ == "__main__":
    main(out=parse_out(__doc__))
