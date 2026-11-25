"""The time-distance diagram of T-Rex, from its tof model.

:func:`main` writes, at the default control parameters:

- ``trex_time_distance.png``: the time-distance diagram of two pulses, from
  :func:`.time_distance_diagram`;
- ``trex_time_distance_params.txt``: the parameters of that run, from
  :func:`.time_distance.parameters_table`.

It needs only mcdgs and tof. Run the module for them::

    pixi run python -m mcdgs.diagnostics.trex.time_distance --out FOLDER
"""

from pathlib import Path

from mcdgs.diagnostics import cli
from mcdgs.diagnostics.time_distance import parameters_table
from mcdgs.instruments.trex.model import Params, TRex
from mcdgs.tof import time_distance_diagram

N_PULSES = 2
"""Pulses of the diagram."""
VISIBLE_RAYS = 2000
"""Most rays drawn per pulse, of the neutrons that reach the sample."""
BLOCKED_RAYS = 5000


# --8<-- [start:instrument]
def define_instrument() -> TRex:
    """Return the instrument of this diagnostic: :class:`.TRex` at its default parameters."""
    params = Params()  # the default control parameters
    return TRex(params=params)


# --8<-- [end:instrument]
"""Most rays drawn per pulse, in gray, of the neutrons a chopper stops."""


def control_parameters(trex: TRex) -> list[tuple[str, str]]:
    """Return the control parameters of T-Rex that set the choppers, as formatted rows.

    Parameters
    ----------
    trex:
        The instrument.

    Returns
    -------
    :
        ``mode``, ``lambda0``, ``rrm`` and ``t0_offset``, each with its value and unit.
    """
    p = trex.params
    return [
        ("`mode`", f'`"{p.mode}"`'),
        ("`lambda0`", f"{p.lambda0:g} Å"),
        ("`rrm`", f"{p.rrm}"),
        ("`t0_offset`", f"{p.t0_offset * 1e3:g} ms"),
    ]


def main(out: Path | None = None, neutrons: int = 1_000_000, seed: int = 1) -> list[Path]:
    """Draw the time-distance diagram of T-Rex and list its parameters, at the defaults.

    Parameters
    ----------
    out:
        Folder to write to. If ``None``, ``trex_time_distance`` in the results folder of
        :meth:`.Workspace.find`.
    neutrons:
        Number of tof neutrons per pulse.
    seed:
        Seed of the random numbers.

    Returns
    -------
    :
        The paths of the diagram and of the table of its parameters.
    """
    folder = cli.output_folder(out, "trex_time_distance")

    # --8<-- [start:time_distance]
    trex = define_instrument()
    figure = time_distance_diagram(
        trex,
        neutrons,
        seed=seed,
        n_pulses=N_PULSES,
        visible_rays=VISIBLE_RAYS,
        blocked_rays=BLOCKED_RAYS,
    )
    path = cli.save_figure(figure, folder / "trex_time_distance.png")
    # --8<-- [end:time_distance]

    table = parameters_table(
        trex,
        control_parameters(trex),
        neutrons=neutrons,
        seed=seed,
        n_pulses=N_PULSES,
        visible_rays=VISIBLE_RAYS,
        blocked_rays=BLOCKED_RAYS,
    )
    return [path, cli.save_text(table, folder / "trex_time_distance_params.txt")]


if __name__ == "__main__":
    main(out=cli.parse_out(__doc__))
