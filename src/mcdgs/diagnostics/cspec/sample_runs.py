"""What the CSPEC diagnostics with a sample share: the two parts of a run split at the sample.

A diagnostic with a sample, such as :mod:`.diagnostics.cspec.diffraction`, runs McStas in two
parts with :class:`.Simulation`: from the source to ``SAMPLE``, which holds ``MONITOR`` and dumps
the beam there into a beam file, then from that beam through the sample to the detector panels.
This module reduces both parts, with :mod:`mcdgs.reduction`:

1. :func:`timing_monitor` is ``mc_monitor``, the monitor after the monochromatic choppers, with
   time bins of ``MC_MONITOR_TIME_RESOLUTION``, finer than those of CSPEC, to time each
   sub-pulse: the instrument of each diagnostic puts it in place of its own.
2. :func:`subpulse_timing` finds the sub-pulses at ``MONITOR`` of the first part with
   :func:`.find_patches`. CSPEC has no beamstop monitor, so, as in
   :mod:`.diagnostics.cspec.incident_energy`, the incident energy of each comes from its flight
   from the source to ``MONITOR`` alone, leaving the source at ``t0_offset``, and the time it
   reaches the sample from the flight on, with :func:`.time_at_sample`.
   :func:`subpulse_timing_of_run` does it for the folder of a run.
3. :func:`read_detector` reads the neutrons of the twelve detector panels of the second part, by
   pixel. That part runs in the frame of ``SAMPLE``, which is the beam frame: the origin is the
   sample, z along the beam, ``BEAM_DIRECTION``, y up and x to the left of the beam.
4. :func:`pixel_solid_angle` gives the solid angle of each pixel of the flat panels, for
   :func:`.normalized_events`, which gives each neutron its energy and momentum transfer, per
   incident neutron and per steradian.

It is not a diagnostic itself: it has no instrument of its own.
"""

from pathlib import Path

import numpy as np
import scipp as sc

from mcdgs.core.components import Monitor_nD
from mcdgs.diagnostics.cspec.incident_energy import MC_MONITOR_TIME_RESOLUTION
from mcdgs.diagnostics.elastic_resolution import BEAM_DIRECTION
from mcdgs.instruments.cspec.detectors import (
    PANEL_HEIGHT,
    PANEL_PIXELS,
    PANEL_WIDTH,
    PANEL_X_BINS,
    PANEL_Y_BINS,
    detector_panels,
    panel_angle,
)
from mcdgs.instruments.cspec.model import CSPEC, Params
from mcdgs.instruments.cspec.monitors import mc_monitor
from mcdgs.instruments.cspec.sources import wavelength_band
from mcdgs.mcstas import Run, read_events, read_mcstas_maps
from mcdgs.reduction.binning import monitor_bins
from mcdgs.reduction.incident_energy import incident_flight
from mcdgs.reduction.inelastic import time_at_sample
from mcdgs.reduction.patches import MIN_FRACTION, find_patches
from mcdgs.reduction.solid_angle import rectangle_solid_angle

__all__ = [
    "BEAM_DIRECTION",
    "MONITOR",
    "SAMPLE",
    "pixel_solid_angle",
    "read_detector",
    "subpulse_timing",
    "subpulse_timing_of_run",
    "timing_monitor",
]

MONITOR = "mc_monitor"  # times the sub-pulses, before the sample
SAMPLE = "sample_pos"


def timing_monitor(params: Params) -> Monitor_nD:
    """Create ``mc_monitor`` with time bins of ``MC_MONITOR_TIME_RESOLUTION``.

    Parameters
    ----------
    params:
        The control parameters of the instrument: the monitor spans the band of ``lambda0``.

    Returns
    -------
    :
        A new monitor, in the place of ``mc_monitor``.
    """
    return mc_monitor(
        lambda_limits=wavelength_band(params.lambda0), time_resolution=MC_MONITOR_TIME_RESOLUTION
    )


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


def subpulse_timing_of_run(to_sample: str | Path, instrument: CSPEC) -> sc.DataArray:
    """Find the sub-pulses of a run to the sample, without running McStas.

    Parameters
    ----------
    to_sample:
        The folder of the run to the sample.
    instrument:
        The instrument of the run, as it ran, as :func:`subpulse_timing`.

    Returns
    -------
    :
        The sub-pulses of :func:`subpulse_timing`.
    """
    return subpulse_timing(Run.load(Path(to_sample)), instrument)


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
