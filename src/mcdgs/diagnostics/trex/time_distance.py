"""The time-distance diagram of T-Rex, from its tof model.

:func:`main` writes, at the default control parameters:

- ``trex_time_distance.png``: the time-distance diagram of two pulses, from
  :func:`.time_distance_diagram`;
- ``trex_time_distance_params.txt``: the parameters of that run, from :func:`parameters_table`.

It needs only mcdgs and tof. Run the module for them::

    pixi run python -m mcdgs.diagnostics.trex.time_distance --out FOLDER
"""

from pathlib import Path

from mcdgs.core.components import NXdisk_chopper
from mcdgs.core.utils import evaluate
from mcdgs.diagnostics.cli import markdown_table, output_folder, parse_out, save_figure, save_text
from mcdgs.instruments.trex.model import TRex
from mcdgs.mcstas.run_folders import default_output_dir
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
        (
            "Source",
            f"ESS pulse of beam port {source.sector}{source.beamline} at {source.frequency:g} Hz",
        ),
        ("Band of the source", f"{lmin:.3f} to {lmax:.3f} Å"),
        ("`neutrons`, per pulse", f"{neutrons:,}"),
        ("`n_pulses`", n_pulses),
        ("`seed`", seed),
        ("`visible_rays`, per pulse", f"{visible_rays:,}"),
        ("`blocked_rays`, per pulse", f"{blocked_rays:,}"),
    ]
    p = trex.params
    control = [
        ("`mode`", f'`"{p.mode}"`'),
        ("`lambda0`", f"{p.lambda0:g} Å"),
        ("`rrm`", p.rrm),
        ("`t0_offset`", f"{p.t0_offset * 1e3:g} ms"),
    ]
    choppers = [
        (
            f"`{name}`",
            f"{distance:.3f} m",
            f"{float(evaluate(component.nu, params)):g} Hz",
            f"{float(evaluate(component.delay, params)) * 1e3:.3f} ms",
        )
        for name, component, distance in after_source
        if isinstance(component, NXdisk_chopper)
    ]
    return "\n".join(
        [
            markdown_table(["tof run", "Value"], run),
            markdown_table(["Control parameter", "Value"], control),
            markdown_table(["Chopper", "Distance from the source", "Frequency", "Delay"], choppers),
        ]
    )


def main(out: Path | None = None, neutrons: int = 1_000_000, seed: int = 1) -> list[Path]:
    """Draw the time-distance diagram of T-Rex and list its parameters, at the defaults.

    Parameters
    ----------
    out:
        Folder to write to. If ``None``, ``mcstas_output/trex_time_distance``, beside the McStas
        runs.
    neutrons:
        Number of tof neutrons per pulse.
    seed:
        Seed of the random numbers.

    Returns
    -------
    :
        The paths of the diagram and of the table of its parameters.
    """
    folder = output_folder(out, default_output_dir() / "trex_time_distance")

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
    path = save_figure(figure, folder / "trex_time_distance.png")
    # --8<-- [end:time_distance]

    table = parameters_table(
        trex,
        neutrons=neutrons,
        seed=seed,
        n_pulses=N_PULSES,
        visible_rays=VISIBLE_RAYS,
        blocked_rays=BLOCKED_RAYS,
    )
    return [path, save_text(table, folder / "trex_time_distance_params.txt")]


if __name__ == "__main__":
    main(out=parse_out(__doc__))
