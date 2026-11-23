"""What the T-Rex diagnostics with a sample share: the two parts of a run split at the sample.

A diagnostic with a sample, such as :mod:`.diagnostics.trex.elastic_plane`, runs McStas in two
parts with :class:`.Simulation`: from the source to ``SAMPLE``, which holds the monitors of
``MONITORS`` and dumps the beam there into a beam file, then from that beam through the sample to
``DETECTOR``. This module reduces both parts, with :mod:`mcdgs.reduction`:

1. :func:`subpulse_timing` finds the sub-pulses at the pre-sample and beamstop monitors of the
   first part with :func:`.find_patches`: both come before ``SAMPLE`` in beam order. It gives each
   sub-pulse its incident energy, from its flight between the monitors with
   :func:`.subpulse_flight`, and the time it reaches the sample, with :func:`.time_at_sample`.
   :func:`subpulse_timing_of_run` does it for the folder of a run, on the monitor bins of its
   instrument.
2. :func:`read_detector` reads the neutrons of the banana detector of the second part, by pixel.
   That part runs in the frame of ``SAMPLE``, which is the beam frame: the origin is the sample,
   z the axis of ``arm_mc``, straight along the beam from the pre-sample monitor to the sample,
   ``BEAM_DIRECTION``, y up and x to the left of the beam.
3. :func:`pixel_solid_angle` gives the solid angle of each pixel, for
   :func:`.normalized_events`, which gives each neutron its energy and momentum transfer, per
   incident neutron and per steradian.

It is not a diagnostic itself: it has no instrument of its own.
"""

from pathlib import Path

import numpy as np
import scipp as sc

from mcdgs.instruments.trex.detectors import (
    BANANA_HEIGHT,
    BANANA_RADIUS,
    BANANA_THETA_BINS,
    BANANA_THETA_LIMITS,
    BANANA_Y_BINS,
)
from mcdgs.instruments.trex.model import TRex
from mcdgs.mcstas import Run, read_events, read_mcstas_maps
from mcdgs.reduction.binning import monitor_bins
from mcdgs.reduction.inelastic import time_at_sample
from mcdgs.reduction.patches import MIN_FRACTION, find_patches, subpulse_flight
from mcdgs.reduction.solid_angle import rectangle_solid_angle

__all__ = [
    "BEAM_DIRECTION",
    "DETECTOR",
    "MONITORS",
    "SAMPLE",
    "pixel_solid_angle",
    "read_detector",
    "subpulse_timing",
    "subpulse_timing_of_run",
]

MONITORS = ("pre_sample_monitor", "beamstop_monitor")  # in beam order
DETECTOR = "banana_detector"
SAMPLE = "sample_pos"
BEAM_DIRECTION = (0.0, 0.0, 1.0)  # the incident beam in the frame of sample_pos: along arm_mc


def subpulse_timing(
    run: Run,
    bins: dict[str, dict[str, sc.Variable]],
    positions: dict[str, list[float]],
    min_fraction: float = MIN_FRACTION,
) -> sc.DataArray:
    """Find the sub-pulses at both monitors, and their incident energy and time at the sample.

    The sub-pulses are the patches of :func:`.find_patches` with at least ``min_fraction`` of
    the intensity of the strongest, paired by :func:`.subpulse_flight`: the n-th at one monitor
    is the n-th at the other.

    Parameters
    ----------
    run:
        The McStas run to the sample, with the monitors of ``MONITORS``.
    bins:
        The bin edges of the time and wavelength map of each monitor of ``MONITORS``, by name,
        as :func:`.monitor_bins`.
    positions:
        The position of every component, in m, by name, as :meth:`.Instrument.positions`.
    min_fraction:
        The weakest patch that counts as a sub-pulse, as a fraction of the strongest.

    Returns
    -------
    :
        The intensity of each sub-pulse at the pre-sample monitor, along ``subpulse``, with the
        coordinates of :func:`.incident_flight` between the monitors (``Ltotal``, ``tof``,
        ``incident_wavelength``, ``incident_energy``), ``time_at_pre_sample_monitor`` and
        ``time_at_sample``.

    Raises
    ------
    ValueError
        If the monitors do not have as many sub-pulses, or one arrives at the beamstop monitor
        first.
    """
    maps = read_mcstas_maps(run, {name: bins[name] for name in MONITORS})
    pre_sample, beamstop = (find_patches(maps[name].values, bins[name]) for name in MONITORS)
    pre_sample_position = sc.vector(positions["pre_sample_monitor"], unit="m")
    pulses = subpulse_flight(
        pre_sample,
        pre_sample_position,
        beamstop,
        sc.vector(positions["beamstop_monitor"], unit="m"),
        min_fraction=min_fraction,
    )
    for name in ("time_at_end", "wavelength_at_start"):
        del pulses.coords[name]
    pulses.coords["time_at_pre_sample_monitor"] = pulses.coords.pop("time_at_start")
    pulses.coords["time_at_sample"] = time_at_sample(
        pulses.coords["time_at_pre_sample_monitor"],
        pre_sample_position,
        sc.vector(positions[SAMPLE], unit="m"),
        pulses.coords["incident_energy"],
    )
    return pulses


def subpulse_timing_of_run(to_sample: str | Path, instrument: TRex) -> sc.DataArray:
    """Find the sub-pulses of a run to the sample, without running McStas.

    Parameters
    ----------
    to_sample:
        The folder of the run to the sample.
    instrument:
        The instrument of the run, as it ran: it gives the bins of the monitors of ``MONITORS``
        and the positions of the components.

    Returns
    -------
    :
        The sub-pulses of :func:`subpulse_timing`.
    """
    bins = {monitor: monitor_bins(getattr(instrument, monitor).bin_limits) for monitor in MONITORS}
    return subpulse_timing(Run.load(Path(to_sample)), bins, instrument.positions())


def read_detector(from_sample: str | Path, entry: int = 1) -> sc.DataArray:
    """Read the neutrons of the banana detector of a run from the sample, by pixel.

    The run has no source: the incident beam is given by its direction, ``BEAM_DIRECTION``.

    Parameters
    ----------
    from_sample:
        The folder of the run from the sample.
    entry:
        The entry to read, from 1: of a scan, the scan point.

    Returns
    -------
    :
        The neutrons of :func:`.read_events`, binned by pixel, in n/s.
    """
    return read_events(
        from_sample, DETECTOR, "source", SAMPLE, beam_direction=BEAM_DIRECTION, entry=entry
    )


def pixel_solid_angle(binned: sc.DataArray) -> sc.Variable:
    """Return the solid angle of each pixel of the banana detector, seen from the sample.

    A pixel is ``BANANA_THETA_BINS`` of the arc of ``BANANA_THETA_LIMITS`` wide, at
    ``BANANA_RADIUS``, and ``BANANA_Y_BINS`` of ``BANANA_HEIGHT`` high, on a cylinder around the
    vertical through the sample, as :func:`.rectangle_solid_angle`.

    Parameters
    ----------
    binned:
        The neutrons by pixel, as :func:`read_detector`, with the pixel ``position`` and the
        ``sample_position``.

    Returns
    -------
    :
        The solid angle of each pixel, in sr, along the pixel dimension.
    """
    theta_start, theta_stop = BANANA_THETA_LIMITS
    width = BANANA_RADIUS * np.radians(theta_stop - theta_start) / BANANA_THETA_BINS
    return rectangle_solid_angle(
        binned.coords["position"],
        binned.coords["sample_position"],
        sc.scalar(width, unit="m"),
        sc.scalar(BANANA_HEIGHT / BANANA_Y_BINS, unit="m"),
        shape="cylinder",
    )
