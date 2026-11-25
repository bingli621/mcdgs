"""The time-distance diagram of CSPEC, from its tof model.

tof takes the sample monitor of CSPEC as a detector, so the rays reach the sample.

:func:`main` writes, at the default control parameters:

- ``cspec_time_distance.png``: the time-distance diagram of two pulses, from
  :func:`.time_distance_diagram`;
- ``cspec_time_distance_params.txt``: the parameters of that run, from
  :func:`.time_distance.parameters_table`.

It needs only mcdgs and tof. Run the module for them::

    pixi run python -m mcdgs.diagnostics.cspec.time_distance --out FOLDER
"""

from pathlib import Path

from mcdgs.diagnostics import cli
from mcdgs.diagnostics.time_distance import parameters_table
from mcdgs.instruments.cspec.model import CSPEC, Params
from mcdgs.tof import time_distance_diagram

N_PULSES = 2
"""Pulses of the diagram."""
VISIBLE_RAYS = 2000
"""Most rays drawn per pulse, of the neutrons that reach the sample."""
BLOCKED_RAYS = 5000
"""Most rays drawn per pulse, in gray, of the neutrons a chopper stops."""


# --8<-- [start:instrument]
def define_instrument() -> CSPEC:
    """Return the instrument of this diagnostic: :class:`.CSPEC` at its default parameters."""
    params = Params()  # the default control parameters
    return CSPEC(params=params)


# --8<-- [end:instrument]


def control_parameters(cspec: CSPEC) -> list[tuple[str, str]]:
    """Return the control parameters of CSPEC that set the choppers, as formatted rows.

    Parameters
    ----------
    cspec:
        The instrument.

    Returns
    -------
    :
        ``lambda0``, ``nu_mc``, ``n``, ``t0_offset`` and ``bw3_delay``, each with its value and
        unit.
    """
    p = cspec.params
    return [
        ("`lambda0`", f"{p.lambda0:g} Å"),
        ("`nu_mc`", f"{p.nu_mc:g} Hz"),
        ("`n`", f"{p.n}"),
        ("`t0_offset`", f"{p.t0_offset * 1e3:g} ms"),
        ("`bw3_delay`", f"{p.bw3_delay * 1e3:g} ms"),
    ]


def main(out: Path | None = None, neutrons: int = 1_000_000, seed: int = 1) -> list[Path]:
    """Draw the time-distance diagram of CSPEC and list its parameters, at the defaults.

    Parameters
    ----------
    out:
        Folder to write to. If ``None``, ``cspec_time_distance`` in the results folder of
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
    folder = cli.output_folder(out, "cspec_time_distance")

    # --8<-- [start:time_distance]
    cspec = define_instrument()
    figure = time_distance_diagram(
        cspec,
        neutrons,
        seed=seed,
        n_pulses=N_PULSES,
        visible_rays=VISIBLE_RAYS,
        blocked_rays=BLOCKED_RAYS,
    )
    path = cli.save_figure(figure, folder / "cspec_time_distance.png")
    # --8<-- [end:time_distance]

    table = parameters_table(
        cspec,
        control_parameters(cspec),
        neutrons=neutrons,
        seed=seed,
        n_pulses=N_PULSES,
        visible_rays=VISIBLE_RAYS,
        blocked_rays=BLOCKED_RAYS,
    )
    return [path, cli.save_text(table, folder / "cspec_time_distance_params.txt")]


if __name__ == "__main__":
    main(out=cli.parse_out(__doc__))
