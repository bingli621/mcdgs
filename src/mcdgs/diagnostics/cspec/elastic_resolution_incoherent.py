"""The elastic resolution of CSPEC, from vanadium: the intensity of each sub-pulse against E.

:class:`CSPECElasticResolution` is :class:`.CSPEC` with the vanadium cylinder of :func:`.vanadium`
on ``sample_pos``, and the detector panels. Its control parameters are those of :class:`.CSPEC`,
so the choppers turn and the beam on the sample is a train of sub-pulses, each with its own
incident energy. Vanadium scatters incoherently and elastically, in every direction, so every
neutron it scatters has an energy transfer of zero: the width of its elastic line is the energy
resolution of the instrument, and whatever lies away from it comes from the reduction or the
instrument, not the sample. The cylinder scatters only once, into the window of
:func:`detector_window`, which covers the panels.

:func:`main` runs it in two parts, with :class:`.Simulation`, as
:mod:`.diagnostics.trex.elastic_resolution_incoherent` does: from the source to ``sample_pos``,
which dumps the beam there into a beam file, then from that beam through the vanadium to the
detector, into the folders of :func:`run_names`.

The reduction is that of T-Rex, with the sub-pulses of :mod:`.diagnostics.cspec.incident_energy`:

1. :func:`subpulse_timing` finds the sub-pulses at ``mc_monitor``, the monitor after the
   monochromatic choppers, of the first part. CSPEC has no beamstop monitor, so their incident
   energy comes from their flight from the source to ``mc_monitor`` alone, every sub-pulse
   leaving the source at ``t0_offset``, and the time they reach the sample from the flight on
   from ``mc_monitor``.
2. :func:`read_detector` reads the neutrons of the twelve detector panels of the second part, by
   pixel, and :func:`pixel_solid_angle` gives the solid angle of each pixel of the flat panels.
3. :func:`.normalized_events` keeps the neutrons inside the arrival window of a sub-pulse, gives
   them their energy transfer, and divides each by the solid angle of its pixel and by the flux
   of its sub-pulse.
4. :func:`.energy_spectra` sums them into a histogram of the energy transfer of each sub-pulse,
   averaged over the solid angle of the pixels: per incident neutron, per steradian and per meV.

:func:`.fit_elastic_lines` fits a Gaussian to the elastic line of each sub-pulse, against the
final energy Ef, and gives its FWHM. :func:`.plot_spectra` draws the spectra with their errors and
their Gaussians, against E, and :func:`.plot_q_e_maps` the S(Q, E) map of each sub-pulse.

:func:`main` is :func:`simulate`, which runs McStas, then :func:`analyze`, which reduces and
draws. Run it to do both::

    pixi run python -c "from mcdgs.diagnostics.cspec.elastic_resolution_incoherent import main; main()"

Running the module only reduces the runs again, without McStas, with :func:`analyze`::

    pixi run python -m mcdgs.diagnostics.cspec.elastic_resolution_incoherent
"""

import math
from pathlib import Path
from typing import NamedTuple

import numpy as np
import scipp as sc

from mcdgs.core.components import Incoherent
from mcdgs.diagnostics import cli
from mcdgs.diagnostics.cspec.incident_energy import MC_MONITOR_TIME_RESOLUTION
from mcdgs.diagnostics.elastic_resolution import (
    BEAM_DIRECTION,
    energy_spectra,
    fit_elastic_lines,
    lines_table,
    plot_q_e_maps,
    plot_spectra,
    q_e_maps,
    two_theta,
)
from mcdgs.instruments.cspec.detectors import (
    PANEL_ARM_ANGLES,
    PANEL_DISTANCE,
    PANEL_HEIGHT,
    PANEL_PIXELS,
    PANEL_WIDTH,
    PANEL_X_BINS,
    PANEL_Y_BINS,
    detector_panels,
    panel_angle,
)
from mcdgs.instruments.cspec.model import CSPEC, DEFAULT_PARAMS, Params
from mcdgs.instruments.cspec.monitors import mc_monitor
from mcdgs.instruments.cspec.sources import wavelength_band
from mcdgs.mcstas import Run, Simulation, Workspace, read_events, read_mcstas_maps
from mcdgs.reduction.binning import monitor_bins
from mcdgs.reduction.incident_energy import incident_flight
from mcdgs.reduction.inelastic import normalized_events, time_at_sample
from mcdgs.reduction.patches import MIN_FRACTION, find_patches
from mcdgs.reduction.solid_angle import rectangle_solid_angle
from mcdgs.samples import vanadium

NAME = "cspec_elastic_resolution_incoherent"
MONITOR = "mc_monitor"  # times the sub-pulses, before the sample
SAMPLE = "sample_pos"
FOCUS_MARGIN = 1.0  # deg, the focusing window of the sample is this much wider than the detector


def detector_window() -> dict[str, float]:
    """Return the focusing window that covers the twelve detector panels, seen from the sample.

    In the horizontal plane, the window runs from the outer edge of the panel at the smallest
    angle, on one side of the beam, to that of the panel at the largest, on the other. It is
    centered in the horizontal plane, which makes McStas span it exactly in height, and reaches
    the top and the bottom edges of the panels at their centers, where they are nearest and so
    the highest in angle. It is ``FOCUS_MARGIN`` wider on every side: the sample is not a point.

    Returns
    -------
    :
        ``target_x``, ``target_z``, ``focus_aw`` and ``focus_ah`` of the sample, in the frame
        of ``sample_pos``.

    Examples
    --------
    >>> window = detector_window()
    >>> round(window["focus_aw"], 2), round(window["focus_ah"], 2)
    (177.28, 55.13)
    """
    half_width = math.degrees(math.atan2(PANEL_WIDTH / 2, PANEL_DISTANCE))
    angles = [panel_angle(i) for i in range(len(PANEL_ARM_ANGLES))]
    start, stop = min(angles) - half_width, max(angles) + half_width
    elevation = math.degrees(math.atan2(PANEL_HEIGHT / 2, PANEL_DISTANCE))
    azimuth = math.radians((start + stop) / 2)
    return {
        "target_x": math.sin(azimuth),
        "target_z": math.cos(azimuth),
        "focus_aw": stop - start + 2 * FOCUS_MARGIN,
        "focus_ah": 2 * (elevation + FOCUS_MARGIN),
    }


# --8<-- [start:instrument]
def sample() -> Incoherent:
    """Return the vanadium cylinder of :func:`.vanadium`, scattering once into the detector.

    It scatters into the window of :func:`detector_window` only.
    """
    return vanadium(name="sample", order=1, **detector_window())


class CSPECElasticResolution(CSPEC):
    """CSPEC with a vanadium cylinder, from the source to the detector panels.

    The cylinder is that of :func:`sample`. ``mc_monitor`` has time bins of
    ``MC_MONITOR_TIME_RESOLUTION``, finer than those of CSPEC, to time each sub-pulse.

    Parameters
    ----------
    params:
        Control parameters.
    source_masked:
        As in :class:`.CSPEC`. True by default: the masked source of
        :func:`.sources.masked_source` spends the rays where the choppers let neutrons through.
    """

    name = "cspec_elastic_resolution"  # the McStas instrument name

    def __init__(self, params: Params = DEFAULT_PARAMS, source_masked: bool = True) -> None:
        super().__init__(sample=sample(), params=params, source_masked=source_masked)
        self.mc_monitor = mc_monitor(  # in its place, with finer time bins
            lambda_limits=wavelength_band(params.lambda0),
            time_resolution=MC_MONITOR_TIME_RESOLUTION,
        )


def define_instrument() -> CSPECElasticResolution:
    """Return the instrument of this diagnostic: :class:`CSPECElasticResolution` at its defaults."""
    params = Params()  # the default control parameters
    return CSPECElasticResolution(params=params)


# --8<-- [end:instrument]


def subpulse_timing(
    run: Run, instrument: CSPEC, min_fraction: float = MIN_FRACTION
) -> sc.DataArray:
    """Find the sub-pulses at ``mc_monitor``, and their incident energy and time at the sample.

    The sub-pulses are the patches of :func:`.find_patches` with at least ``min_fraction`` of the
    intensity of the strongest. The incident energy comes from the flight from the source to the
    monitor, through :func:`.incident_flight`, every sub-pulse leaving the source at
    ``t0_offset``, and the time at the sample from the flight on, with :func:`.time_at_sample`.

    Parameters
    ----------
    run:
        The McStas run to the sample, with ``mc_monitor``.
    instrument:
        The instrument of the run, as it ran: it gives the bins of ``mc_monitor``, the positions
        of the components and ``t0_offset``.
    min_fraction:
        The weakest patch that counts as a sub-pulse, as a fraction of the strongest.

    Returns
    -------
    :
        The intensity of each sub-pulse at ``mc_monitor``, in n/s, along ``subpulse``, with the
        coordinates of :func:`.incident_flight` from the source (``Ltotal``, ``tof``,
        ``incident_wavelength``, ``incident_energy``), ``time_at_mc_monitor`` and
        ``time_at_sample``.
    """
    bins = {MONITOR: monitor_bins(getattr(instrument, MONITOR).bin_limits)}
    intensity = read_mcstas_maps(run, bins)[MONITOR]
    subpulses = find_patches(intensity.values, bins[MONITOR]).strongest(min_fraction)
    positions = instrument.positions()
    monitor_position = sc.vector(positions[MONITOR], unit="m")
    at_monitor = sc.array(dims=["subpulse"], values=subpulses.time, unit="ms")
    flight = incident_flight(
        sc.vector(positions["source"], unit="m"),
        sc.scalar(instrument.params.t0_offset, unit="s").to(unit="ms"),
        monitor_position,
        at_monitor,
    )
    pulses = sc.DataArray(
        sc.array(dims=["subpulse"], values=subpulses.intensity, unit="1/s"),  # n/s
        coords={"time_at_mc_monitor": at_monitor, **flight},
    )
    pulses.coords["time_at_sample"] = time_at_sample(
        at_monitor,
        monitor_position,
        sc.vector(positions[SAMPLE], unit="m"),
        pulses.coords["incident_energy"],
    )
    return pulses


def read_detector(from_sample: str | Path) -> sc.DataArray:
    """Read the neutrons of the twelve detector panels of a run from the sample, by pixel.

    The run has no source: the incident beam is given by its direction, ``BEAM_DIRECTION``.

    Parameters
    ----------
    from_sample:
        The folder of the run from the sample.

    Returns
    -------
    :
        The neutrons of :func:`.read_events` of every panel, binned by pixel, in n/s, the panels
        one after the other along the pixel dimension.
    """
    panels = [
        read_events(from_sample, name, "source", SAMPLE, beam_direction=BEAM_DIRECTION)
        for name in detector_panels()
    ]
    (dim,) = panels[0].dims
    return sc.concat(panels, dim)


def pixel_solid_angle(binned: sc.DataArray) -> sc.Variable:
    """Return the solid angle of each pixel of the detector panels, seen from the sample.

    A pixel is ``PANEL_X_BINS`` of ``PANEL_WIDTH`` wide and ``PANEL_Y_BINS`` of ``PANEL_HEIGHT``
    high, on a flat panel: its normal is that of its panel, at :func:`.panel_angle` from the beam,
    found from its number, ``PANEL_PIXELS`` per panel, as :func:`.rectangle_solid_angle`.

    Parameters
    ----------
    binned:
        The neutrons by pixel, as :func:`read_detector`, with the pixel ``position``,
        ``pixel_id`` and the ``sample_position``.

    Returns
    -------
    :
        The solid angle of each pixel, in sr, along the pixel dimension.
    """
    panel = binned.coords["pixel_id"].values // PANEL_PIXELS
    angle = np.radians([panel_angle(int(i)) for i in panel])
    (dim,) = binned.dims
    normal = sc.vectors(
        dims=[dim], values=np.stack([np.sin(angle), np.zeros_like(angle), np.cos(angle)], axis=1)
    )
    return rectangle_solid_angle(
        binned.coords["position"],
        binned.coords["sample_position"],
        sc.scalar(PANEL_WIDTH / PANEL_X_BINS, unit="m"),
        sc.scalar(PANEL_HEIGHT / PANEL_Y_BINS, unit="m"),
        normal=normal,
    )


class Reduction(NamedTuple):
    """The reduction of a run of :func:`simulate`, from :func:`reduce_runs`."""

    pulses: sc.DataArray  # the sub-pulses of :func:`subpulse_timing`
    spectra: sc.DataArray  # of :func:`.energy_spectra`
    maps: sc.DataArray  # of :func:`.q_e_maps`


def run_names() -> tuple[str, str]:
    """Return the run folders: the part to ``sample_pos``, and the part from it.

    Examples
    --------
    >>> run_names()[0]
    'cspec_elastic_resolution_incoherent_run_to_sample'
    """
    return f"{NAME}_run_to_sample", f"{NAME}_run_from_sample"


def simulate(ncount: float = 1e9, seed: int = 1) -> tuple[Path, Path]:
    """Run :class:`CSPECElasticResolution` in McStas, in two parts split at ``sample_pos``.

    The first part runs from the source to ``sample_pos``, with ``mc_monitor`` and the beam file;
    the second from that beam through the vanadium to the detector panels, both into the folders
    of :func:`run_names`, in the runs folder of :meth:`.Workspace.find`. Each replaces the folder
    of an earlier run.

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
    to_name, from_name = run_names()
    sim = Simulation(define_instrument())
    to_sample = sim.run_to_sample(folder_name=to_name, ncount=ncount, seed=seed)
    # sample_pos -> detector, from the beam file, in the frame of sample_pos: the beam frame
    from_sample = sim.run_from_sample(folder_name=from_name, seed=seed + 1)
    return to_sample.folder, from_sample.folder


def reduce_runs(
    to_sample: str | Path | None = None, from_sample: str | Path | None = None
) -> Reduction:
    """Reduce the two parts of a McStas run of :func:`simulate`, without running McStas.

    Every pixel is kept: the direct beam passes between the two panels closest to it.

    Parameters
    ----------
    to_sample, from_sample:
        The folders of the first and the second part. If ``None``, those of :func:`run_names`
        in the runs folder of :meth:`.Workspace.find`, where :func:`simulate` writes.

    Returns
    -------
    :
        The sub-pulses of :func:`subpulse_timing`, the spectra of :func:`.energy_spectra` and
        the maps of :func:`.q_e_maps`.
    """
    workspace = Workspace.find()
    to_name, from_name = run_names()
    to_sample = Path(to_sample or workspace.run_folder(to_name))
    from_sample = Path(from_sample or workspace.run_folder(from_name))
    pulses = subpulse_timing(Run.load(to_sample), define_instrument())  # as simulate runs it
    binned = read_detector(from_sample)
    pixels = pixel_solid_angle(binned)
    events = normalized_events(binned, pulses, pixels)
    return Reduction(
        pulses,
        energy_spectra(events, pulses, pixels.sum()),
        q_e_maps(events, pulses, two_theta(binned), pixels),
    )


def analyze(out: Path | None = None) -> list[Path]:
    """Reduce the runs, and save their spectra, their fits and their S(Q, E) maps.

    It runs no McStas, so it can run again on the same runs; :func:`main` runs them first. The
    files are ``cspec_elastic_resolution_incoherent.png``, from :func:`.plot_spectra`,
    ``cspec_elastic_resolution_incoherent_lines.txt``, the table of :func:`.lines_table`, and
    ``cspec_elastic_resolution_incoherent_sqe.png``, from :func:`.plot_q_e_maps`.

    Parameters
    ----------
    out:
        Folder of the files. If ``None``, ``cspec_elastic_resolution_incoherent`` in the results
        folder of :meth:`.Workspace.find`.

    Returns
    -------
    :
        The paths of the files.
    """
    folder = cli.output_folder(out, NAME)
    pulses, spectra, maps = reduce_runs()
    conditions = define_instrument().running_conditions()
    lines = fit_elastic_lines(spectra)
    table = lines_table(lines)
    print(table, end="")
    files = [cli.save_text(table, folder / f"{NAME}_lines.txt")]
    fig = plot_spectra(
        spectra,
        lines,
        f"CSPEC, vanadium: the energy transfer of {pulses.sizes['subpulse']} sub-pulses\n"
        f"{conditions}",
    )
    files.append(cli.save_figure(fig, folder / f"{NAME}.png"))
    fig = plot_q_e_maps(maps, f"CSPEC, vanadium: S(Q, E) of each sub-pulse\n{conditions}")
    files.append(cli.save_figure(fig, folder / f"{NAME}_sqe.png"))
    return files


def main(out: Path | None = None, ncount: float = 1e9, seed: int = 1) -> list[Path]:
    """Run McStas with :func:`simulate`, then reduce and draw with :func:`analyze`.

    Parameters
    ----------
    out:
        Folder of the files, as :func:`analyze`.
    ncount, seed:
        As :func:`simulate`.

    Returns
    -------
    :
        The paths of the files.
    """
    simulate(ncount=ncount, seed=seed)
    return analyze(out)


if __name__ == "__main__":
    analyze(out=cli.parse_out(__doc__))  # McStas runs only with main
