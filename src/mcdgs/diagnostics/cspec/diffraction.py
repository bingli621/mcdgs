"""Diffraction from a La2Ni7 crystal on CSPEC: the horizontal scattering plane, in (Qx, Qz).

Every neutron in the arrival window of its sub-pulse counts, whatever its energy transfer: the
inelastic scattering is weak next to the Bragg peaks, so the plane is a diffraction pattern.

:class:`CSPECWithLa2Ni7` is :class:`.CSPEC` with the La2Ni7 crystal of :func:`.la2ni7` on
``sample_pos``, as on T-Rex, and the detector panels, which list every neutron reaching them. Its
control parameters are those of :class:`.CSPEC`. Few neutrons reach a 5 mm crystal, so McStas
repeats each one that does ``SAMPLE_SPLIT`` times, with its weight divided among them.

:func:`main` runs it in two parts, with :class:`.Simulation`: from the source to ``sample_pos``,
which dumps the beam there into a beam file, then from that beam through the crystal to the
detector. The second part can run again from the same beam, without the long first part.

The reduction is that of :mod:`.diagnostics.cspec.sample_runs`, with :mod:`mcdgs.reduction`:

1. :func:`.sample_runs.subpulse_timing` finds the sub-pulses at ``mc_monitor`` of the first
   part, with their incident energy, from their flight from the source leaving at
   ``t0_offset``, and the time they reach the sample.
2. :func:`.sample_runs.read_detector` reads the neutrons of the twelve detector panels of the
   second part, by pixel, in the frame of ``sample_pos``, which is the beam frame: the origin is
   the sample, z along the beam, y up and x to the left of the beam.
3. :func:`.events_to_q_and_energy`, within :func:`.normalized_events`, gives each pixel and
   sub-pulse the arrival window of :func:`.arrival_window`, keeps the neutrons that
   :func:`.assign_subpulse` finds inside the window of a sub-pulse, whatever their energy
   transfer within it, and converts them with ``QE_GRAPH`` to the momentum transfer (Qx, Qy,
   Qz).
4. :func:`.normalized_events` divides each neutron by the solid angle of its pixel, from
   :func:`.sample_runs.pixel_solid_angle`, and by the incoming flux of its sub-pulse, its
   intensity at ``mc_monitor``: per incident neutron and per steradian.

:func:`.plot_plane` draws the intensity in (Qx, Qz) of every sub-pulse, summed over a slice of Qy
near the horizontal plane, on a log color scale ``LOG_DECADES`` deep.

:func:`main` is :func:`simulate`, which runs McStas, then :func:`analyze`, which reduces and
draws. Run it to do both::

    pixi run python -c "from mcdgs.diagnostics.cspec.diffraction import main; main()"

Running the module only reduces the runs again, without McStas, with :func:`analyze`;
:func:`reduce_runs` gives the sub-pulses and the neutrons alone::

    pixi run python -m mcdgs.diagnostics.cspec.diffraction
"""

from pathlib import Path

import scipp as sc

from mcdgs.diagnostics import cli
from mcdgs.diagnostics.cspec import sample_runs as runs
from mcdgs.diagnostics.diffraction import plot_plane, print_subpulse_table
from mcdgs.instruments.cspec.model import CSPEC, DEFAULT_PARAMS, Params
from mcdgs.mcstas import Simulation, Workspace
from mcdgs.reduction.inelastic import normalized_events
from mcdgs.samples import la2ni7

NAME = "cspec_diffraction"
TO_SAMPLE = f"{NAME}_run_to_sample"  # the run folder of the first part, source -> sample_pos
FROM_SAMPLE = f"{NAME}_run_from_sample"  # that of the second part, sample_pos -> detector
SAMPLE_SPLIT = 10  # McStas SPLIT on the crystal: each neutron reaching it is repeated this often
# the color scale spans this many decades below the highest bin, one fewer than for T-Rex: the arcs
# of all the sub-pulses meet near the origin, so the highest bin is far above the incoherent arcs
LOG_DECADES = 3


# --8<-- [start:instrument]
class CSPECWithLa2Ni7(CSPEC):
    """CSPEC with a La2Ni7 crystal, from the source to the detector panels.

    The crystal is that of :func:`.la2ni7`, split ``SAMPLE_SPLIT`` times. ``mc_monitor`` is
    :func:`.sample_runs.timing_monitor`, with finer time bins, to time each sub-pulse.

    Parameters
    ----------
    params:
        Control parameters.
    source_masked:
        As in :class:`.CSPEC`. True by default: the masked source of
        :func:`.sources.masked_source` spends the rays where the choppers let neutrons through.
    """

    name = "cspec_la2ni7"  # the McStas instrument name

    def __init__(self, params: Params = DEFAULT_PARAMS, source_masked: bool = True) -> None:
        sample = la2ni7(name="sample", split=SAMPLE_SPLIT)
        super().__init__(sample=sample, params=params, source_masked=source_masked)
        self.mc_monitor = runs.timing_monitor(params)  # in its place, with finer time bins


def define_instrument() -> CSPECWithLa2Ni7:
    """Return the instrument of this diagnostic: :class:`CSPECWithLa2Ni7` at its defaults."""
    params = Params()  # the default control parameters
    return CSPECWithLa2Ni7(params=params)


# --8<-- [end:instrument]


def simulate(ncount: float = 1e9, seed: int = 1) -> tuple[Path, Path]:
    """Run :class:`CSPECWithLa2Ni7` in McStas, in two parts split at ``sample_pos``.

    The first part runs from the source to ``sample_pos`` into ``TO_SAMPLE``, with
    ``mc_monitor`` and the beam file; the second from that beam through the crystal to the
    detector panels, into ``FROM_SAMPLE``, both in the runs folder of :meth:`.Workspace.find`.
    Each replaces the folder of an earlier run.

    Parameters
    ----------
    ncount:
        Number of McStas rays of the first part, from the masked source. The second part runs
        every neutron of the beam file.
    seed:
        Seed of the random numbers of the first part; the second part takes ``seed + 1``.

    Returns
    -------
    :
        The folders of the first and the second part.
    """
    sim = Simulation(define_instrument())
    to_sample = sim.run_to_sample(folder_name=TO_SAMPLE, ncount=ncount, seed=seed)
    # sample_pos -> detector, from the beam file, in the frame of sample_pos: the beam frame
    from_sample = sim.run_from_sample(folder_name=FROM_SAMPLE, seed=seed + 1)
    return to_sample.folder, from_sample.folder


def reduce_runs(
    to_sample: str | Path | None = None, from_sample: str | Path | None = None
) -> tuple[sc.DataArray, sc.DataArray]:
    """Reduce the two parts of a McStas run of :func:`simulate`, without running McStas.

    Every pixel is kept: the direct beam passes between the two panels closest to it.

    Parameters
    ----------
    to_sample, from_sample:
        The folders of the first and the second part. If ``None``, ``TO_SAMPLE`` and
        ``FROM_SAMPLE`` in the runs folder of :meth:`.Workspace.find`, where :func:`simulate`
        writes.

    Returns
    -------
    :
        The sub-pulses of :func:`.sample_runs.subpulse_timing`, and the neutrons of
        :func:`.normalized_events`.
    """
    workspace = Workspace.find()
    to_sample = Path(to_sample or workspace.run_folder(TO_SAMPLE))
    from_sample = Path(from_sample or workspace.run_folder(FROM_SAMPLE))
    # the first part: the incident energy and the time at the sample of each sub-pulse
    pulses = runs.subpulse_timing_of_run(to_sample, define_instrument())  # as simulate runs it
    binned = runs.read_detector(from_sample)
    return pulses, normalized_events(binned, pulses, runs.pixel_solid_angle(binned))


def analyze(
    out: Path | None = None,
    to_sample: str | Path | None = None,
    from_sample: str | Path | None = None,
) -> Path:
    """Reduce the two parts of a McStas run, print the sub-pulses and save the (Qx, Qz) plane.

    It runs no McStas, so it can run again on the same runs; :func:`main` runs them first. The
    figure is ``cspec_diffraction.png``, from :func:`.plot_plane`.

    Parameters
    ----------
    out:
        Folder of the figure. If ``None``, ``cspec_diffraction`` in the results folder of
        :meth:`.Workspace.find`.
    to_sample, from_sample:
        The folders of the two parts, as :func:`reduce_runs`.

    Returns
    -------
    :
        The path of the figure.
    """
    pulses, events = reduce_runs(to_sample, from_sample)
    print_subpulse_table(pulses, events)

    fig = plot_plane(
        events,
        f"CSPEC, La2Ni7: (Qx, Qz) of {pulses.sizes['subpulse']} sub-pulses\n"
        f"{define_instrument().running_conditions()}",
        log_decades=LOG_DECADES,
    )
    return cli.save_figure(fig, cli.output_folder(out, NAME) / f"{NAME}.png")


def main(out: Path | None = None, ncount: float = 1e9, seed: int = 1) -> Path:
    """Run McStas with :func:`simulate`, then reduce and draw it with :func:`analyze`.

    Parameters
    ----------
    out:
        Folder of the figure, as :func:`analyze`.
    ncount, seed:
        As :func:`simulate`.

    Returns
    -------
    :
        The path of the figure.
    """
    to_sample, from_sample = simulate(ncount=ncount, seed=seed)
    return analyze(out, to_sample, from_sample)


if __name__ == "__main__":
    analyze(out=cli.parse_out(__doc__))  # McStas runs only with main
