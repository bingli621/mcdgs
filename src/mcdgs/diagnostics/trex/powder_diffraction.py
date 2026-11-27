"""Diffraction from an Al2O3 powder on T-Rex: its Debye-Scherrer rings, in (Qx, Qz) and in |Q|.

Every neutron in the arrival window of its sub-pulse counts, whatever its energy transfer: the
powder of :func:`.al2o3` scatters elastically only, on its Bragg cones and incoherently.

:class:`TRexWithAl2O3` is :class:`.TRex` with the Al2O3 powder of :func:`.al2o3` on
``sample_pos``, with the beamstop and the banana detector, which lists every neutron reaching
it. Its control parameters are those of :class:`.TRex`. The powder spreads each Bragg cone
within ``d_phi`` of the horizontal plane only, the height of the window of
:func:`.solid_angles.detector_window`, which covers the detector. McStas repeats each neutron
reaching the powder ``SAMPLE_SPLIT`` times, with its weight divided among them. ``slit_mc`` is
closed to the size of the powder, as wide as the cylinder and as high, so the pre-sample monitor
counts the beam that reaches the powder, not the wider beam of the guide.

:func:`simulate` runs it in two parts, with :class:`.Simulation`: from the source to
``sample_pos``, which dumps the beam there into a beam file, then from that beam through the
powder to the detector. The beam before the sample goes through the closed slit, so it is not
that of :mod:`.diagnostics.trex.diffraction`.

The reduction is that of :mod:`.diagnostics.trex.diffraction`, with the sub-pulses of the first
part, and the same neutrons and normalization: per incident neutron and per steradian.

:func:`.plot_plane` draws the intensity in (Qx, Qz) of every sub-pulse, summed over a slice of
Qy near the horizontal plane, ``|Qy| < QY_HALF_WIDTH``, on a log color scale. Each sub-pulse has
its own incident wavevector ki, so its elastic neutrons lie near its own circle, of radius ki
through the origin, and its Bragg neutrons where that circle crosses the rings of the powder,
of radius 2π / d about the origin.

A powder scatters the same way about the beam, so its pattern depends on |Q| alone:
:func:`.powder_patterns` reduces the neutrons of each sub-pulse to its intensity in the elastic
|Q| of their pixel, from every pixel, over the solid angle of the pixels that see each bin, and
:func:`.plot_powder_patterns` draws those of all the sub-pulses over each other.

:func:`.detector_views` puts the same neutrons on the detector, in the horizontal angle θ and the
height of their pixel, and :func:`.plot_detector_views` draws the view of each sub-pulse on its
own panel: each Bragg cone of the powder crosses the detector as a Debye-Scherrer ring.

:func:`main` is :func:`simulate`, which runs McStas, then :func:`analyze`, which reduces and
draws. Run it to do both::

    pixi run python -c "from mcdgs.diagnostics.trex.powder_diffraction import main; main()"

Running the module only reduces the runs again, without McStas, with :func:`analyze`;
:func:`reduce_runs` gives the sub-pulses, the neutrons, the patterns
and the detector views alone::

    pixi run python -m mcdgs.diagnostics.trex.powder_diffraction
"""

from pathlib import Path

import scipp as sc

from mcdgs.diagnostics import cli
from mcdgs.diagnostics.diffraction import plot_plane, print_subpulse_table
from mcdgs.diagnostics.elastic_resolution import two_theta
from mcdgs.diagnostics.powder_diffraction import (
    detector_views,
    plot_detector_views,
    plot_powder_patterns,
    powder_patterns,
)
from mcdgs.diagnostics.trex import sample_runs as runs
from mcdgs.diagnostics.trex.solid_angles import detector_window
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.mcstas import Simulation, Workspace
from mcdgs.reduction.inelastic import normalized_events
from mcdgs.samples import al2o3

NAME = "trex_powder_diffraction"
TO_SAMPLE = f"{NAME}_run_to_sample"  # the run folder of the first part, source -> sample_pos
FROM_SAMPLE = f"{NAME}_run_from_sample"  # that of the second part, sample_pos -> detector
SAMPLE_SPLIT = 10  # McStas SPLIT on the powder: each neutron reaching it is repeated this often
QY_HALF_WIDTH = 0.25  # 1/Å: the (Qx, Qz) plane keeps the neutrons with |Qy| below this
VIEW_HEADROOM = 5.0  # the top of the color scale of the detector views, over the highest bin


# --8<-- [start:instrument]
class TRexWithAl2O3(TRex):
    """T-Rex with an Al2O3 powder, from the source to the banana detector.

    The powder is that of :func:`.al2o3`, its Bragg cones within ``d_phi`` of the horizontal
    plane, the height of the window of :func:`.solid_angles.detector_window`, split
    ``SAMPLE_SPLIT`` times. ``slit_mc`` is closed to the size of the powder: as wide as the
    cylinder, ``2 * radius``, and as high, ``yheight``.

    Parameters
    ----------
    params:
        Control parameters.
    source_masked:
        As in :class:`.TRex`. True by default: the masked source of
        :func:`.sources.masked_source` spends the rays where the choppers let neutrons through.
    """

    name = "trex_al2o3"  # the McStas instrument name

    def __init__(self, params: Params = DEFAULT_PARAMS, source_masked: bool = True) -> None:
        sample = al2o3(name="sample", d_phi=detector_window()["focus_ah"], split=SAMPLE_SPLIT)
        super().__init__(sample=sample, params=params, source_masked=source_masked)
        self.slit_mc.xwidth = 2 * sample.radius  # the slit, closed to the size of the powder
        self.slit_mc.yheight = sample.yheight


def define_instrument() -> TRexWithAl2O3:
    """Return the instrument of this diagnostic: :class:`TRexWithAl2O3` at its defaults."""
    params = Params()  # the default control parameters
    return TRexWithAl2O3(params=params)


# --8<-- [end:instrument]


def simulate(ncount: float = 1e9, seed: int = 1) -> tuple[Path, Path]:
    """Run :class:`TRexWithAl2O3` in McStas, in two parts split at ``sample_pos``.

    The first part runs from the source to ``sample_pos`` into ``TO_SAMPLE``, with the monitors
    of the incident energy and the beam file; the second from that beam through the powder to
    the detector, into ``FROM_SAMPLE``, both in the runs folder of :meth:`.Workspace.find`. Each
    replaces the folder of an earlier run.

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
) -> tuple[sc.DataArray, sc.DataArray, sc.DataArray, sc.DataArray]:
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
        The sub-pulses of :func:`.sample_runs.subpulse_timing`, the neutrons of
        :func:`.normalized_events`, the patterns of :func:`.powder_patterns`, from every
        pixel, and the detector views of :func:`.detector_views`.
    """
    workspace = Workspace.find()
    to_sample = Path(to_sample or workspace.run_folder(TO_SAMPLE))
    from_sample = Path(from_sample or workspace.run_folder(FROM_SAMPLE))
    # the first part: the incident energy and the time at the sample of each sub-pulse
    pulses = runs.subpulse_timing_of_run(to_sample, define_instrument())  # as simulate runs it
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
    figures are ``trex_powder_diffraction.png``, the (Qx, Qz) plane of :func:`.plot_plane`, and
    ``trex_powder_diffraction_q.png``, the intensity in |Q| of every sub-pulse of
    :func:`.plot_powder_patterns`, and ``trex_powder_diffraction_detector.png``, the detector
    view of every sub-pulse of :func:`.plot_detector_views`.

    Parameters
    ----------
    out:
        Folder of the figures. If ``None``, ``trex_powder_diffraction`` in the results folder of
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
        f"T-Rex, Al2O3 powder: (Qx, Qz) of {n} sub-pulses\n{conditions}",
        qy_half_width=QY_HALF_WIDTH,
    )
    pattern = plot_powder_patterns(
        patterns,
        f"T-Rex, Al2O3 powder: I(|Q|) of each of {n} sub-pulses\n{conditions}",
        monitor="pre-sample monitor",
    )
    detector = plot_detector_views(
        views,
        f"T-Rex, Al2O3 powder: detector view of each of {n} sub-pulses\n{conditions}",
        monitor="pre-sample monitor",
        headroom=VIEW_HEADROOM,
    )
    return [
        cli.save_figure(plane, folder / f"{NAME}.png"),
        cli.save_figure(pattern, folder / f"{NAME}_q.png"),
        cli.save_figure(detector, folder / f"{NAME}_detector.png"),
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
    analyze(out=cli.parse_out(__doc__))  # McStas runs only with main
