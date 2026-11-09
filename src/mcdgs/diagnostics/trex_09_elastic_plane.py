"""The horizontal scattering plane of T-Rex with a La2Ni7 crystal, in (Qx, Qz).

:class:`TRexWithLa2Ni7` is :class:`.TRex` with the La2Ni7 crystal of :func:`.la2ni7` on
``sample_pos``, with the beamstop and the banana detector, which lists every neutron reaching
it. Its control parameters are those of :class:`.TRex`. Few neutrons reach a 5 mm crystal, so
McStas repeats each one that does ``SAMPLE_SPLIT`` times, with its weight divided among them.

:func:`main` runs it in two parts, with :class:`.Simulation`: from the source to ``sample_pos``,
which dumps the beam there into a beam file, then from that beam through the crystal to the
detector. The second part can run again from the same beam, without the long first part.

The reduction follows :mod:`mcdgs.reduction`:

1. :func:`subpulse_timing` finds the sub-pulses at the pre-sample and beamstop monitors of the
   first part with :func:`.find_patches`, as :mod:`.trex_08_determine_ei` does: both come before
   ``sample_pos`` in beam order. It gives each sub-pulse its incident energy, from its flight
   between the monitors with :func:`.subpulse_flight`, and the time it reaches the sample, with
   :func:`.time_at_sample`.
2. :func:`.read_events` reads the neutrons of the banana detector of the second part, by pixel.
   That part runs in the frame of ``sample_pos``, which is the beam frame: the origin is the
   sample, z the axis of ``arm_mc``, straight along the beam from the pre-sample monitor to the
   sample, y up and x to the left of the beam.
3. :func:`.events_to_q_and_energy`, within :func:`normalized_events`, gives each pixel and
   sub-pulse the arrival window of :func:`.arrival_window`, keeps the neutrons that
   :func:`.assign_subpulse` finds inside the window of a sub-pulse, whatever their energy
   transfer within it, and converts them with ``QE_GRAPH`` to the momentum transfer (Qx, Qy,
   Qz).
4. :func:`normalized_events` divides each neutron by the solid angle of its pixel, from
   :func:`pixel_solid_angle`, and by the incoming flux of its sub-pulse, its intensity at the
   pre-sample monitor: per incident neutron and per steradian.

:func:`plot_plane` draws the intensity in (Qx, Qz) of every sub-pulse, summed over a slice of
Qy near the horizontal plane, ``|Qy| < QY_HALF_WIDTH``, on a log color scale. Each sub-pulse has
its own incident wavevector ki, so its elastic neutrons lie near its own circle, of radius ki
through the origin, where the detector covers it.
Run the module to do it::

    pixi run python -m mcdgs.diagnostics.trex_09_elastic_plane

:func:`main` is :func:`simulate`, which runs McStas, then :func:`analyze`, which reduces and
draws. Reduce the runs again without McStas with :func:`analyze`, or :func:`reduce_runs` for the
sub-pulses and the neutrons alone::

    pixi run python -c "from mcdgs.diagnostics.trex_09_elastic_plane import analyze; analyze()"
"""

from pathlib import Path

import matplotlib.pyplot as plt
import mcstasscript as ms
import numpy as np
import scipp as sc
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure
from mcstasscript.data.data import McStasData

from mcdgs.diagnostics.cli import output_folder
from mcdgs.instruments.trex.detectors import (
    BANANA_HEIGHT,
    BANANA_RADIUS,
    BANANA_THETA_BINS,
    BANANA_THETA_LIMITS,
    BANANA_Y_BINS,
)
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.mcstas.nexus_io import read_events
from mcdgs.mcstas.run_folders import default_output_dir
from mcdgs.mcstas.simulation import Simulation
from mcdgs.reduction.binning import bin_edges
from mcdgs.reduction.inelastic import events_to_q_and_energy, time_at_sample
from mcdgs.reduction.patches import MIN_FRACTION, find_patches, subpulse_flight
from mcdgs.reduction.solid_angle import rectangle_solid_angle
from mcdgs.samples import la2ni7
from mcdgs.tof.reader import monitor_bins

MONITORS = ("pre_sample_monitor", "beamstop_monitor")  # in beam order
DETECTOR = "banana_detector"
SAMPLE = "sample_pos"
NAME = "trex_elastic_plane"
TO_SAMPLE = f"{NAME}_run_to_sample"  # the run folder of the first part, source -> sample_pos
FROM_SAMPLE = f"{NAME}_run_from_sample"  # that of the second part, sample_pos -> detector
BEAM_DIRECTION = (0.0, 0.0, 1.0)  # the incident beam in the frame of sample_pos: along arm_mc
Q_STEP = 0.02  # 1/Å, the bin width in Qx and Qz
LOG_DECADES = 4  # the color scale spans this many decades below the highest bin
QY_HALF_WIDTH = 0.5  # 1/Å: the plane keeps the neutrons with |Qy| below this, near horizontal
DIRECT_BEAM_RADIUS = 0.2  # 1/Å: the direct beam, far above the rest, is left out of the scale
SAMPLE_SPLIT = 100  # McStas SPLIT on the crystal: each neutron reaching it is repeated this often


class TRexWithLa2Ni7(TRex):
    """T-Rex with a La2Ni7 crystal, from the source to the banana detector.

    The crystal is that of :func:`.la2ni7`, split ``SAMPLE_SPLIT`` times.

    Parameters
    ----------
    params:
        Control parameters.
    """

    name = "trex_la2ni7"  # the McStas instrument name

    def __init__(self, params: Params = DEFAULT_PARAMS) -> None:
        sample = la2ni7(name="sample", split=SAMPLE_SPLIT)
        super().__init__(sample=sample, params=params)


def subpulse_timing(
    mcstas_data: list[McStasData],
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
    mcstas_data:
        The monitors of the McStas run, as ``mcstasscript.load_data`` reads them.
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
    pre_sample, beamstop = (
        # McStasScript has the map as (wavelength, time), the y axis first
        find_patches(
            np.asarray(ms.name_search(name, mcstas_data).Intensity, dtype=float).T, bins[name]
        )
        for name in MONITORS
    )
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


def pixel_solid_angle(binned: sc.DataArray) -> sc.Variable:
    """Return the solid angle of each pixel of the banana detector, seen from the sample.

    A pixel is ``BANANA_THETA_BINS`` of the arc of ``BANANA_THETA_LIMITS`` wide, at
    ``BANANA_RADIUS``, and ``BANANA_Y_BINS`` of ``BANANA_HEIGHT`` high, on a cylinder around the
    vertical through the sample, as :func:`.rectangle_solid_angle`.

    Parameters
    ----------
    binned:
        The neutrons by pixel, as :func:`.read_events`, with the pixel ``position`` and the
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


def normalized_events(binned: sc.DataArray, pulses: sc.DataArray) -> sc.DataArray:
    """Reduce the neutrons of the detector, per incident neutron and per steradian.

    Each neutron is divided by the solid angle of its pixel, from :func:`pixel_solid_angle`,
    then reduced with :func:`.events_to_q_and_energy`, and divided by the intensity of its
    sub-pulse at the pre-sample monitor, its incoming flux.

    Parameters
    ----------
    binned:
        The neutrons by pixel, as :func:`.read_events`, in n/s.
    pulses:
        The sub-pulses of :func:`subpulse_timing`, with their intensity at the pre-sample
        monitor, in n/s.

    Returns
    -------
    :
        The neutrons of :func:`.events_to_q_and_energy`, in 1/sr: per neutron of their
        sub-pulse at the pre-sample monitor. Each keeps the ``solid_angle`` of its pixel, so
        that times it gives the fraction of the incident neutrons it stands for.
    """
    solid_angle = pixel_solid_angle(binned)
    binned = binned.copy()
    binned.bins.coords["solid_angle"] = sc.bins_like(binned, solid_angle)
    events = events_to_q_and_energy(binned / solid_angle, pulses)
    flux = pulses.data[events.coords["subpulse"].values]  # the flux of each neutron's sub-pulse
    events.data = events.data / sc.array(dims=events.dims, values=flux.values, unit=flux.unit)
    return events


def print_subpulse_table(pulses: sc.DataArray, events: sc.DataArray) -> None:
    """Print the incident energy, ki and time at the sample of each sub-pulse.

    Parameters
    ----------
    pulses:
        The sub-pulses of :func:`subpulse_timing`.
    events:
        The neutrons of :func:`normalized_events`.
    """
    print(
        "  sub-pulse  t_sample [ms]  λi [Å]  Ei [meV]  ki [1/Å]  flux [n/s]  neutrons"
        "  on the detector"
    )
    subpulse = events.coords["subpulse"].values
    n = pulses.sizes["subpulse"]
    counts = np.bincount(subpulse, minlength=n)
    # the fraction of the incident neutrons that reach the detector, in a window
    fractions = np.bincount(
        subpulse, events.values * events.coords["solid_angle"].values, minlength=n
    )
    for i, (t, wavelength, ei, flux, count, fraction) in enumerate(
        zip(
            pulses.coords["time_at_sample"].to(unit="ms").values,
            pulses.coords["incident_wavelength"].to(unit="Å").values,
            pulses.coords["incident_energy"].to(unit="meV").values,
            pulses.values,
            counts,
            fractions,
        )
    ):
        k = 2 * np.pi / wavelength
        print(
            f"  {i:9d}  {t:13.3f}  {wavelength:6.4f}  {ei:8.3f}  {k:8.4f}  {flux:10.4g}"
            f"  {count:8d}  {fraction:15.4g}"
        )


def plot_plane(
    events: sc.DataArray,
    title: str,
    step: float = Q_STEP,
    qy_half_width: float = QY_HALF_WIDTH,
) -> Figure:
    """Draw the intensity in (Qx, Qz) of every sub-pulse, near the horizontal plane.

    Only the neutrons with ``|Qy| < qy_half_width`` count, summed over Qy: the detector also
    sees reflections well above and below the plane, which a sum over all Qy would project
    onto it.

    It is drawn on the log color scale of :func:`draw_plane`.

    Parameters
    ----------
    events:
        The neutrons of :func:`.events_to_q_and_energy`.
    title:
        The title of the figure.
    step:
        The bin width in Qx and Qz, in 1/Å.
    qy_half_width:
        Half the width of the slice in Qy, in 1/Å.

    Returns
    -------
    :
        The figure.
    """
    qy = events.coords["Qy"].to(unit="1/Å").values
    (dim,) = events.dims
    events = events[dim, np.flatnonzero(np.abs(qy) < qy_half_width)]
    edges = {}
    for name in ("Qx", "Qz"):
        values = events.coords[name].to(unit="1/Å").values
        start, stop = np.floor(values.min() * 10) / 10, np.ceil(values.max() * 10) / 10
        edges[name] = bin_edges(name, start, stop, step, unit="1/Å")
    plane = events.hist(Qz=edges["Qz"], Qx=edges["Qx"])  # summed over the slice in Qy
    return draw_plane(plane, f"{title}\n|Qy| < {qy_half_width:g} 1/Å")


def beyond_direct_beam_max(
    plane: sc.DataArray, direct_beam_radius: float = DIRECT_BEAM_RADIUS
) -> float:
    """Return the highest bin of a (Qz, Qx) histogram beyond the direct beam.

    Parameters
    ----------
    plane:
        The histogram, with the bin edges ``Qx`` and ``Qz``, in 1/Å, and any dimensions before
        them.
    direct_beam_radius:
        The |Q| of the direct beam, in 1/Å, left out.

    Returns
    -------
    :
        The highest value of the bins with ``|Q| > direct_beam_radius``.
    """
    qx, qz = (sc.midpoints(plane.coords[name]).to(unit="1/Å").values for name in ("Qx", "Qz"))
    outside = np.hypot(*np.meshgrid(qz, qx, indexing="ij")) > direct_beam_radius
    return float(plane.values[..., outside].max())


def draw_plane(
    plane: sc.DataArray,
    title: str,
    xlabel: str = "Qx [1/Å]",
    ylabel: str = "Qz [1/Å], along the beam",
    direct_beam_radius: float = DIRECT_BEAM_RADIUS,
    vmax: float | None = None,
) -> Figure:
    """Draw a histogram in (Qx, Qz) on a log color scale.

    The scale is turbo, ``LOG_DECADES`` deep below the highest bin beyond the direct beam,
    ``|Q| > direct_beam_radius``, with the empty bins in white: the direct beam that passes the
    beamstop, a few bins far above all others, would leave the rest at the bottom of the scale,
    so it takes the color of the top instead. The axes have the same scale, so that circles look
    round.

    Parameters
    ----------
    plane:
        The histogram, with the bin edges ``Qx`` and ``Qz``, in 1/Å, evenly spaced.
    title:
        The title of the figure.
    xlabel, ylabel:
        The labels of the axes.
    direct_beam_radius:
        The |Q| of the direct beam, in 1/Å, that the color scale leaves out.
    vmax:
        The top of the color scale. If ``None``, the highest bin beyond the direct beam.

    Returns
    -------
    :
        The figure.
    """
    step = float(np.diff(plane.coords["Qx"].values)[0])
    if vmax is None:
        vmax = beyond_direct_beam_max(plane, direct_beam_radius)

    fig, ax = plt.subplots(figsize=(9, 8), layout="constrained")
    fig.suptitle(title)
    mesh = ax.pcolormesh(
        plane.coords["Qx"].values,
        plane.coords["Qz"].values,
        np.ma.masked_less_equal(plane.values, 0),  # empty bins in white
        cmap=plt.get_cmap("turbo").with_extremes(bad="white"),
        norm=LogNorm(vmin=vmax / 10**LOG_DECADES, vmax=vmax, clip=True),
    )
    ax.set(xlabel=xlabel, ylabel=ylabel, aspect="equal")
    ax.grid(True, alpha=0.3)
    fig.colorbar(
        mesh,
        ax=ax,
        shrink=0.8,
        label=f"Intensity [{plane.unit} per ({step:g} 1/Å)²]",
        extend="max",
    )
    return fig


def simulate(ncount: float = 1e9, seed: int = 1) -> tuple[Path, Path]:
    """Run :class:`TRexWithLa2Ni7` in McStas, in two parts split at ``sample_pos``.

    The first part runs from the source to ``sample_pos`` into ``TO_SAMPLE``, with the monitors
    of the incident energy and the beam file; the second from that beam through the crystal to
    the detector, into ``FROM_SAMPLE``, both in :func:`.default_output_dir`. Each replaces the
    folder of an earlier run.

    Parameters
    ----------
    ncount:
        Number of McStas neutrons of the first part. The second part runs every neutron of the
        beam file.
    seed:
        Seed of the random numbers of the first part; the second part takes ``seed + 1``.

    Returns
    -------
    :
        The folders of the first and the second part.
    """
    simulation = Simulation(TRexWithLa2Ni7())  # the default control parameters
    beam_file = simulation.run_to_sample(TO_SAMPLE, ncount=ncount, seed=seed)
    # sample_pos -> detector, from the beam file, in the frame of sample_pos: the beam frame
    detector_file = simulation.run_from_sample(FROM_SAMPLE, beam_file=beam_file, seed=seed + 1)
    return beam_file.parent, detector_file.parent


def subpulse_timing_of_run(to_sample: str | Path | None = None) -> sc.DataArray:
    """Find the sub-pulses of the run to the sample of :func:`simulate`, without running McStas.

    Parameters
    ----------
    to_sample:
        The folder of the first part. If ``None``, ``TO_SAMPLE`` in :func:`.default_output_dir`.

    Returns
    -------
    :
        The sub-pulses of :func:`subpulse_timing`.
    """
    to_sample = Path(to_sample or default_output_dir() / TO_SAMPLE)
    instrument = TRexWithLa2Ni7()  # the default control parameters, as simulate runs
    bins = {monitor: monitor_bins(getattr(instrument, monitor)) for monitor in MONITORS}
    return subpulse_timing(ms.load_data(str(to_sample)), bins, instrument.positions())


def reduce_runs(
    to_sample: str | Path | None = None, from_sample: str | Path | None = None
) -> tuple[sc.DataArray, sc.DataArray]:
    """Reduce the two parts of a McStas run of :func:`simulate`, without running McStas.

    Parameters
    ----------
    to_sample, from_sample:
        The folders of the first and the second part. If ``None``, ``TO_SAMPLE`` and
        ``FROM_SAMPLE`` in :func:`.default_output_dir`, where :func:`simulate` writes.

    Returns
    -------
    :
        The sub-pulses of :func:`subpulse_timing`, and the neutrons of
        :func:`.events_to_q_and_energy`.
    """
    from_sample = Path(from_sample or default_output_dir() / FROM_SAMPLE)
    # the first part: the incident energy and the time at the sample of each sub-pulse
    pulses = subpulse_timing_of_run(to_sample)
    # the second part has no source: the incident beam is given by its direction instead
    binned = read_events(from_sample, DETECTOR, "source", SAMPLE, beam_direction=BEAM_DIRECTION)
    return pulses, normalized_events(binned, pulses)


def analyze(
    out: Path | None = None,
    to_sample: str | Path | None = None,
    from_sample: str | Path | None = None,
) -> Path:
    """Reduce the two parts of a McStas run, print the sub-pulses and save the (Qx, Qz) plane.

    It runs no McStas, so it can run again on the same runs; :func:`main` runs them first. The
    figure is ``trex_elastic_plane.png``, from :func:`plot_plane`.

    Parameters
    ----------
    out:
        Folder of the figure. If ``None``, the folder of the second part, next to its
        ``mccode.h5``.
    to_sample, from_sample:
        The folders of the two parts, as :func:`reduce_runs`.

    Returns
    -------
    :
        The path of the figure.
    """
    from_sample = Path(from_sample or default_output_dir() / FROM_SAMPLE)
    pulses, events = reduce_runs(to_sample, from_sample)
    print_subpulse_table(pulses, events)

    fig = plot_plane(
        events,
        f"T-Rex, La2Ni7: (Qx, Qz) of {pulses.sizes['subpulse']} sub-pulses\n"
        f"{TRexWithLa2Ni7().running_conditions()}",
    )
    figure = output_folder(out, from_sample) / f"{NAME}.png"
    fig.savefig(figure, dpi=150)
    plt.close(fig)
    print(f"Saved {figure}")
    return figure


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
    # main(out=parse_out(__doc__))

    analyze()
