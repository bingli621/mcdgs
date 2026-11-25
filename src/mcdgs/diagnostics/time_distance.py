"""The tables of the parameters of a time-distance diagram, for any instrument.

:func:`parameters_table` lists the tof run of a time-distance diagram, the control parameters
that set the choppers, and each chopper, with its distance from the source, frequency and delay.
The time-distance diagnostic of each instrument gives its own control parameters.
"""

from collections.abc import Sequence

from mcdgs.core import Instrument, evaluate
from mcdgs.core.components import NXdisk_chopper
from mcdgs.diagnostics import cli
from mcdgs.tof import beamline


def parameters_table(
    instrument: Instrument,
    control: Sequence[tuple[str, object]],
    *,
    neutrons: int,
    seed: int | None,
    n_pulses: int,
    visible_rays: int,
    blocked_rays: int,
) -> str:
    """Return the parameters of a tof run of an instrument as Markdown tables.

    Three tables: the tof run, with the band of the source; the control parameters that set the
    choppers; and each chopper, with its distance from the source, frequency and delay.

    Parameters
    ----------
    instrument:
        The instrument, at the parameters of the run.
    control:
        The control parameters that set the choppers, as (name, value) rows, already formatted.
    neutrons:
        Number of tof neutrons per pulse.
    seed:
        Seed of the random numbers.
    n_pulses:
        Number of pulses.
    visible_rays, blocked_rays:
        Most rays drawn per pulse, of the neutrons that reach the last detector and that are
        stopped.

    Returns
    -------
    :
        The tables, separated by blank lines.
    """
    source, after_source = beamline(instrument)
    params = instrument.numeric_instr_params()
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
            cli.markdown_table(["tof run", "Value"], run),
            cli.markdown_table(["Control parameter", "Value"], control),
            cli.markdown_table(
                ["Chopper", "Distance from the source", "Frequency", "Delay"], choppers
            ),
        ]
    )
