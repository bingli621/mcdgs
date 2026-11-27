"""Diffraction from a La2Ni7 crystal on T-Rex: the horizontal scattering plane, in (Qx, Qz).

Every neutron in the arrival window of its sub-pulse counts, whatever its energy transfer: the
inelastic scattering is weak next to the Bragg peaks, so the plane is a diffraction pattern.

:class:`TRexWithLa2Ni7` is :class:`.TRex` with the La2Ni7 crystal of :func:`.la2ni7` on
``sample_pos``, with the beamstop and the banana detector, which lists every neutron reaching
it. Its control parameters are those of :class:`.TRex`. Few neutrons reach a 5 mm crystal, so
McStas repeats each one that does ``SAMPLE_SPLIT`` times, with its weight divided among them.

:func:`main` runs it in two parts, with :class:`.Simulation`: from the source to ``sample_pos``,
which dumps the beam there into a beam file, then from that beam through the crystal to the
detector. The second part can run again from the same beam, without the long first part.

The reduction is that of :mod:`.diagnostics.trex.sample_runs`, with :mod:`mcdgs.reduction`:

1. :func:`.sample_runs.subpulse_timing` finds the sub-pulses at the pre-sample and beamstop
   monitors of the first part, as :mod:`.diagnostics.trex.incident_energy` does, with their
   incident energy and the time they reach the sample.
2. :func:`.sample_runs.read_detector` reads the neutrons of the banana detector of the second
   part, by pixel, in the frame of ``sample_pos``, which is the beam frame: the origin is the
   sample, z the axis of ``arm_mc``, straight along the beam from the pre-sample monitor to the
   sample, y up and x to the left of the beam.
3. :func:`.events_to_q_and_energy`, within :func:`.normalized_events`, gives each pixel and
   sub-pulse the arrival window of :func:`.arrival_window`, keeps the neutrons that
   :func:`.assign_subpulse` finds inside the window of a sub-pulse, whatever their energy
   transfer within it, and converts them with ``QE_GRAPH`` to the momentum transfer (Qx, Qy,
   Qz).
4. :func:`.normalized_events` divides each neutron by the solid angle of its pixel, from
   :func:`.sample_runs.pixel_solid_angle`, and by the incoming flux of its sub-pulse, its
   intensity at the pre-sample monitor: per incident neutron and per steradian.

:func:`.plot_plane` draws the intensity in (Qx, Qz) of every sub-pulse, summed over a slice of
Qy near the horizontal plane, ``|Qy| < QY_HALF_WIDTH``, on a log color scale. Each sub-pulse has
its own incident wavevector ki, so its elastic neutrons lie near its own circle, of radius ki
through the origin, where the detector covers it. :func:`.plot_vertical_planes` draws the
intensity in (Qy, Qz) of each sub-pulse on its own panel, summed over ``|Qx| < QX_HALF_WIDTH``.

:func:`main` is :func:`simulate`, which runs McStas, then :func:`analyze`, which reduces and
draws. Run it to do both::

    pixi run python -c "from mcdgs.diagnostics.trex.diffraction import main; main()"

Running the module only reduces the runs again, without McStas, with :func:`analyze`;
:func:`reduce_runs` gives the sub-pulses and the neutrons alone::

    pixi run python -m mcdgs.diagnostics.trex.diffraction
"""

from pathlib import Path

import scipp as sc

from mcdgs.diagnostics import cli
from mcdgs.diagnostics.diffraction import plot_plane, plot_vertical_planes, print_subpulse_table
from mcdgs.diagnostics.trex import sample_runs as runs
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.mcstas import Simulation, Workspace
from mcdgs.reduction.inelastic import normalized_events
from mcdgs.samples import la2ni7

NAME = "trex_diffraction"
TO_SAMPLE = f"{NAME}_run_to_sample"  # the run folder of the first part, source -> sample_pos
FROM_SAMPLE = f"{NAME}_run_from_sample"  # that of the second part, sample_pos -> detector
SAMPLE_SPLIT = 10  # McStas SPLIT on the crystal: each neutron reaching it is repeated this often


# --8<-- [start:instrument]
class TRexWithLa2Ni7(TRex):
    """T-Rex with a La2Ni7 crystal, from the source to the banana detector.

    The crystal is that of :func:`.la2ni7`, split ``SAMPLE_SPLIT`` times.

    Parameters
    ----------
    params:
        Control parameters.
    source_masked:
        As in :class:`.TRex`. True by default: the masked source of
        :func:`.sources.masked_source` spends the rays where the choppers let neutrons through.
    """

    name = "trex_la2ni7"  # the McStas instrument name

    def __init__(self, params: Params = DEFAULT_PARAMS, source_masked: bool = True) -> None:
        sample = la2ni7(name="sample", split=SAMPLE_SPLIT)
        super().__init__(sample=sample, params=params, source_masked=source_masked)


def define_instrument() -> TRexWithLa2Ni7:
    """Return the instrument of this diagnostic: :class:`TRexWithLa2Ni7` at its defaults."""
    params = Params()  # the default control parameters
    return TRexWithLa2Ni7(params=params)


# --8<-- [end:instrument]


def simulate(ncount: float = 1e9, seed: int = 1) -> tuple[Path, Path]:
    """Run :class:`TRexWithLa2Ni7` in McStas, in two parts split at ``sample_pos``.

    The first part runs from the source to ``sample_pos`` into ``TO_SAMPLE``, with the monitors
    of the incident energy and the beam file; the second from that beam through the crystal to
    the detector, into ``FROM_SAMPLE``, both in the runs folder of :meth:`.Workspace.find`. Each
    replaces the folder of an earlier run.

    Parameters
    ----------
    ncount:
        Number of McStas rays of the first part. With the masked source nearly every ray reaches
        the sample: 1e8 rays leave about as many neutrons in the beam file as 1e10 of the plain
        source, and 1e9 about ten times as many. The second part, and every angle of
        :mod:`.diagnostics.trex.diffraction_scan`, runs every neutron of the beam file.
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
) -> list[Path]:
    """Reduce the two parts of a McStas run, print the sub-pulses and save their Q maps.

    It runs no McStas, so it can run again on the same runs; :func:`main` runs them first. The
    figures are ``trex_diffraction.png``, the (Qx, Qz) plane of :func:`.plot_plane`, and
    ``trex_diffraction_qy_qz.png``, the (Qy, Qz) map of each sub-pulse of
    :func:`.plot_vertical_planes`.

    Parameters
    ----------
    out:
        Folder of the figures. If ``None``, ``trex_diffraction`` in the results folder of
        :meth:`.Workspace.find`.
    to_sample, from_sample:
        The folders of the two parts, as :func:`reduce_runs`.

    Returns
    -------
    :
        The paths of the (Qx, Qz) plane and of the (Qy, Qz) maps.
    """
    from_sample = Path(from_sample or Workspace.find().run_folder(FROM_SAMPLE))
    pulses, events = reduce_runs(to_sample, from_sample)
    print_subpulse_table(pulses, events)

    conditions = define_instrument().running_conditions()
    n = pulses.sizes["subpulse"]
    folder = cli.output_folder(out, NAME)
    plane = plot_plane(events, f"T-Rex, La2Ni7: (Qx, Qz) of {n} sub-pulses\n{conditions}")
    maps = plot_vertical_planes(
        events, pulses, f"T-Rex, La2Ni7: (Qy, Qz) of each of {n} sub-pulses\n{conditions}"
    )
    return [
        cli.save_figure(plane, folder / f"{NAME}.png"),
        cli.save_figure(maps, folder / f"{NAME}_qy_qz.png"),
    ]


def main(out: Path | None = None, ncount: float = 1e9, seed: int = 1) -> list[Path]:
    """Run McStas with :func:`simulate`, then reduce and draw it with :func:`analyze`.

    Parameters
    ----------
    out:
        Folder of the figures, as :func:`analyze`.
    ncount, seed:
        As :func:`simulate`.

    Returns
    -------
    :
        The paths of the figures, as :func:`analyze`.
    """
    to_sample, from_sample = simulate(ncount=ncount, seed=seed)
    return analyze(out, to_sample, from_sample)


if __name__ == "__main__":
    analyze(out=cli.parse_out(__doc__))  # McStas runs only with main, about two minutes
