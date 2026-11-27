"""Diffraction from an Al2O3 powder on CSPEC: its Debye-Scherrer rings, in (Qx, Qz) and in |Q|.

Every neutron in the arrival window of its sub-pulse counts, whatever its energy transfer: the
powder of :func:`.al2o3` scatters elastically only, on its Bragg cones and incoherently.

:class:`CSPECWithAl2O3` is :class:`.CSPEC` with the Al2O3 powder of :func:`.al2o3` on
``sample_pos``, as on T-Rex, and the detector panels, which list every neutron reaching them. Its
control parameters are those of :class:`.CSPEC`. The powder spreads each Bragg cone within
``d_phi`` of the horizontal plane only, the height of the window of
:func:`.elastic_resolution_incoherent_short.detector_window`, which covers the panels. McStas
repeats each neutron reaching the powder ``SAMPLE_SPLIT`` times, with its weight divided among
them.

:func:`simulate` runs only the part from ``sample_pos`` to the detector panels, with
:class:`.Simulation`, from the beam file of :mod:`.diagnostics.cspec.diffraction`: the beam before
the sample does not depend on the sample, so the powder takes that of the crystal.

The reduction is that of :mod:`.diagnostics.cspec.diffraction`, with the sub-pulses of its run to
the sample, and the same neutrons and normalization: per incident neutron and per steradian.

:func:`.plot_plane` draws the intensity in (Qx, Qz) of every sub-pulse, summed over a slice of
Qy near the horizontal plane, ``|Qy| < QY_HALF_WIDTH``, on a log color scale ``LOG_DECADES``
deep, as for the crystal. Each sub-pulse has its own incident wavevector ki, so its elastic
neutrons lie near its own circle, of radius ki through the origin, and its Bragg neutrons where
that circle crosses the rings of the powder, of radius 2π / d about the origin.

A powder scatters the same way about the beam, so its pattern depends on |Q| alone:
:func:`.powder_patterns` reduces the neutrons of each sub-pulse to its intensity in the elastic
|Q| of their pixel, from every pixel, over the solid angle of the pixels that see each bin, and
:func:`.plot_powder_patterns` draws those of all the sub-pulses over each other.

:func:`.detector_views` puts the same neutrons on the detector, in the horizontal angle θ and the
height of their pixel, and :func:`.plot_detector_views` draws the view of each sub-pulse on its
own panel: each Bragg cone of the powder crosses the detector as a Debye-Scherrer ring.

:func:`main` is :func:`simulate`, which runs McStas, then :func:`analyze`, which reduces and
draws. Run it to do both::

    pixi run python -c "from mcdgs.diagnostics.cspec.powder_diffraction import main; main()"

Running the module only reduces the runs again, without McStas, with :func:`analyze`;
:func:`reduce_runs` gives the sub-pulses, the neutrons, the patterns
and the detector views alone::

    pixi run python -m mcdgs.diagnostics.cspec.powder_diffraction
"""

from pathlib import Path

import scipp as sc

from mcdgs.diagnostics import cli
from mcdgs.diagnostics.cspec import diffraction as crystal
from mcdgs.diagnostics.cspec import sample_runs as runs
from mcdgs.diagnostics.cspec.diffraction import LOG_DECADES
from mcdgs.diagnostics.cspec.elastic_resolution_incoherent_short import detector_window
from mcdgs.diagnostics.diffraction import plot_plane, print_subpulse_table
from mcdgs.diagnostics.elastic_resolution import two_theta
from mcdgs.diagnostics.powder_diffraction import (
    detector_views,
    plot_detector_views,
    plot_powder_patterns,
    powder_patterns,
)
from mcdgs.instruments.cspec.model import CSPEC, DEFAULT_PARAMS, Params
from mcdgs.mcstas import Simulation, Workspace
from mcdgs.reduction.inelastic import normalized_events
from mcdgs.samples import al2o3

NAME = "cspec_powder_diffraction"
FROM_SAMPLE = f"{NAME}_run_from_sample"  # the run folder, sample_pos -> detector
SAMPLE_SPLIT = 10  # McStas SPLIT on the powder: each neutron reaching it is repeated this often
QY_HALF_WIDTH = 0.25  # 1/Å: the (Qx, Qz) plane keeps the neutrons with |Qy| below this


# --8<-- [start:instrument]
class CSPECWithAl2O3(CSPEC):
    """CSPEC with an Al2O3 powder, from the source to the detector panels.

    The powder is that of :func:`.al2o3`, its Bragg cones within ``d_phi`` of the horizontal
    plane, the height of the window of :func:`.elastic_resolution_incoherent_short.detector_window`,
    split ``SAMPLE_SPLIT`` times. ``mc_monitor`` is :func:`.sample_runs.timing_monitor`, with
    finer time bins, to time each sub-pulse.

    Parameters
    ----------
    params:
        Control parameters.
    source_masked:
        As in :class:`.CSPEC`. True by default: the masked source of
        :func:`.sources.masked_source` spends the rays where the choppers let neutrons through.
    """

    name = "cspec_al2o3"  # the McStas instrument name

    def __init__(self, params: Params = DEFAULT_PARAMS, source_masked: bool = True) -> None:
        sample = al2o3(name="sample", d_phi=detector_window()["focus_ah"], split=SAMPLE_SPLIT)
        super().__init__(sample=sample, params=params, source_masked=source_masked)
        self.mc_monitor = runs.timing_monitor(params)  # in its place, with finer time bins


def define_instrument() -> CSPECWithAl2O3:
    """Return the instrument of this diagnostic: :class:`CSPECWithAl2O3` at its defaults."""
    params = Params()  # the default control parameters
    return CSPECWithAl2O3(params=params)


# --8<-- [end:instrument]


def simulate(beam_run: str = crystal.TO_SAMPLE, seed: int = 2) -> Path:
    """Run :class:`CSPECWithAl2O3` in McStas from ``sample_pos`` to the panels, from a beam file.

    The run goes into ``FROM_SAMPLE`` in the runs folder of :meth:`.Workspace.find`, replacing
    the folder of an earlier run. It runs every neutron of the beam file.

    Parameters
    ----------
    beam_run:
        The name of the run to the sample whose beam to start from: by default that of
        :mod:`.diagnostics.cspec.diffraction`, which :func:`.diffraction.simulate` makes.
    seed:
        Seed of the random numbers: by default that of the second part of
        :func:`.diffraction.simulate`.

    Returns
    -------
    :
        The run folder.
    """
    sim = Simulation(define_instrument())
    # sample_pos -> detector, from the beam file, in the frame of sample_pos: the beam frame
    return sim.run_from_sample(folder_name=FROM_SAMPLE, beam_run=beam_run, seed=seed).folder


def reduce_runs(
    to_sample: str | Path | None = None, from_sample: str | Path | None = None
) -> tuple[sc.DataArray, sc.DataArray, sc.DataArray, sc.DataArray]:
    """Reduce the two parts of a McStas run of :func:`simulate`, without running McStas.

    Every pixel is kept: the direct beam passes between the two panels closest to it.

    Parameters
    ----------
    to_sample, from_sample:
        The run to the sample and the run of :func:`simulate`. If ``None``, the
        ``TO_SAMPLE`` of :mod:`.diagnostics.cspec.diffraction` and ``FROM_SAMPLE``, in the runs
        folder of :meth:`.Workspace.find`.

    Returns
    -------
    :
        The sub-pulses of :func:`.sample_runs.subpulse_timing`, the neutrons of
        :func:`.normalized_events`, the patterns of :func:`.powder_patterns`, from every
        pixel, and the detector views of :func:`.detector_views`.
    """
    workspace = Workspace.find()
    to_sample = Path(to_sample or workspace.run_folder(crystal.TO_SAMPLE))
    from_sample = Path(from_sample or workspace.run_folder(FROM_SAMPLE))
    # the beam is that of diffraction, which its instrument ran: the incident energy and the
    # time at the sample of each sub-pulse
    pulses = runs.subpulse_timing_of_run(to_sample, crystal.define_instrument())
    binned = runs.read_detector(from_sample)
    pixels, angles = runs.pixel_solid_angle(binned), two_theta(binned)
    binned.bins.coords["two_theta"] = sc.bins_like(binned, angles)  # each neutron, its pixel's
    events = normalized_events(binned, pulses, pixels)
    patterns = powder_patterns(events, pulses, angles, pixels)
    views = detector_views(
        events, pulses, binned.coords["position"], binned.coords["sample_position"], pixels
    )
    return pulses, events, patterns, views


def analyze(
    out: Path | None = None,
    to_sample: str | Path | None = None,
    from_sample: str | Path | None = None,
) -> list[Path]:
    """Reduce the two parts of a McStas run, print the sub-pulses and save their figures.

    It runs no McStas, so it can run again on the same runs; :func:`main` runs them first. The
    figures are ``cspec_powder_diffraction.png``, the (Qx, Qz) plane of :func:`.plot_plane`, and
    ``cspec_powder_diffraction_q.png``, the intensity in |Q| of every sub-pulse of
    :func:`.plot_powder_patterns`, and ``cspec_powder_diffraction_detector.png``, the detector
    view of every sub-pulse of :func:`.plot_detector_views`.

    Parameters
    ----------
    out:
        Folder of the figures. If ``None``, ``cspec_powder_diffraction`` in the results folder of
        :meth:`.Workspace.find`.
    to_sample, from_sample:
        As :func:`reduce_runs`.

    Returns
    -------
    :
        The paths of the (Qx, Qz) plane, of the patterns in |Q| and of the detector views.
    """
    pulses, events, patterns, views = reduce_runs(to_sample, from_sample)
    print_subpulse_table(pulses, events)

    conditions = define_instrument().running_conditions()
    n = pulses.sizes["subpulse"]
    folder = cli.output_folder(out, NAME)
    plane = plot_plane(
        events,
        f"CSPEC, Al2O3 powder: (Qx, Qz) of {n} sub-pulses\n{conditions}",
        qy_half_width=QY_HALF_WIDTH,
        log_decades=LOG_DECADES,
    )
    pattern = plot_powder_patterns(
        patterns,
        f"CSPEC, Al2O3 powder: I(|Q|) of each of {n} sub-pulses\n{conditions}",
        monitor="mc_monitor",
    )
    detector = plot_detector_views(
        views,
        f"CSPEC, Al2O3 powder: detector view of each of {n} sub-pulses\n{conditions}",
        monitor="mc_monitor",
    )
    return [
        cli.save_figure(plane, folder / f"{NAME}.png"),
        cli.save_figure(pattern, folder / f"{NAME}_q.png"),
        cli.save_figure(detector, folder / f"{NAME}_detector.png"),
    ]


def main(out: Path | None = None, beam_run: str = crystal.TO_SAMPLE, seed: int = 2) -> list[Path]:
    """Run McStas with :func:`simulate`, then reduce and draw it with :func:`analyze`.

    The beam file of ``beam_run`` must be there: :func:`.diffraction.simulate` makes it.

    Parameters
    ----------
    out:
        Folder of the figures, as :func:`analyze`.
    beam_run, seed:
        As :func:`simulate`.

    Returns
    -------
    :
        The paths of the figures, as :func:`analyze`.
    """
    from_sample = simulate(beam_run=beam_run, seed=seed)
    to_sample = Workspace.find().run_folder(beam_run)
    return analyze(out, to_sample, from_sample)


if __name__ == "__main__":
    analyze(out=cli.parse_out(__doc__))  # McStas runs only with main
