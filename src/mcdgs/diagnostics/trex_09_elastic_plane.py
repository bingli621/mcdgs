"""The horizontal scattering plane of T-Rex with a La2Ni7 crystal, in (Qx, Qz).

:class:`TRexWithLa2Ni7` is :class:`.TRex` with the La2Ni7 crystal of :func:`.la2ni7` on
``sample_pos``, with the beamstop and the banana detector, which lists every neutron reaching
it. Its control parameters are those of :class:`.TRex`. Few neutrons reach a 5 mm crystal, so
McStas repeats each one that does ``SAMPLE_SPLIT`` times, with its weight divided among them.

The reduction follows :mod:`mcdgs.reduction`:

1. :func:`subpulse_timing` finds the sub-pulses at the pre-sample and beamstop monitors of McStas
   with :func:`.find_patches`, as :mod:`.trex_08_determine_ei` does, and gives each its incident
   energy, from its flight between the monitors with :func:`.subpulse_flight`, and the time it
   reaches the sample, with :func:`.time_at_sample`.
2. :func:`read_events_in_beam_frame` reads the neutrons of the banana detector, by pixel, with
   :func:`.read_events`, in the beam frame of :func:`.beam_frame`: z along the beam between the
   pre-sample monitor and the sample, y up.
3. :func:`.events_to_q_and_energy` gives each pixel and sub-pulse the arrival window of
   :func:`.arrival_window`, keeps the neutrons that :func:`.assign_subpulse` finds inside the
   window of a sub-pulse, whatever their energy transfer within it, and converts them with
   ``QE_GRAPH`` to the momentum transfer (Qx, Qy, Qz).

:func:`plot_plane` draws the intensity in (Qx, Qz) of every sub-pulse, summed over Qy and in a
thin slice of Qy. Each sub-pulse has its own incident wavevector ki, so its elastic neutrons lie
on its own circle, of radius ki through the origin, where the detector covers it. Run the module
to do it::

    pixi run python -m mcdgs.diagnostics.trex_09_elastic_plane
"""

from pathlib import Path

import matplotlib.pyplot as plt
import mcstasscript as ms
import numpy as np
import scipp as sc
from matplotlib.axes import Axes
from matplotlib.collections import QuadMesh
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure
from mcstasscript.data.data import McStasData

from mcdgs.diagnostics.cli import output_folder, parse_out
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.mcstas.nexus_io import read_events
from mcdgs.mcstas.simulation import Simulation
from mcdgs.reduction.binning import bin_edges
from mcdgs.reduction.frame import beam_frame
from mcdgs.reduction.inelastic import events_to_q_and_energy, time_at_sample
from mcdgs.reduction.patches import MIN_FRACTION, find_patches, subpulse_flight
from mcdgs.samples import la2ni7
from mcdgs.tof.reader import monitor_bins

MONITORS = ("pre_sample_monitor", "beamstop_monitor")  # in beam order
DETECTOR = "banana_detector"
SAMPLE = "sample_pos"
Q_STEP = 0.02  # 1/Å, the bin width in Qx and Qz
QY_HALF_WIDTH = 0.05  # 1/Å, half the width of the slice in Qy near the horizontal plane
LOG_DECADES = 4  # the log color scale spans this many decades below the highest bin
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


def read_events_in_beam_frame(folder: str | Path) -> sc.DataArray:
    """Read the neutrons of the banana detector of a McStas run, in the beam frame.

    Parameters
    ----------
    folder:
        Run folder with the ``mccode.h5`` to read.

    Returns
    -------
    :
        The neutrons of :func:`.read_events`, with the pre-sample monitor as the source, and
        their positions in the frame of :func:`.beam_frame`.
    """
    # source_position only sets the direction of the incident beam, from it to the sample: the
    # pre-sample monitor, 0.3 m before the sample, gives it where the curved guide has bent it
    binned = read_events(folder, DETECTOR, MONITORS[0], SAMPLE)
    sample_position = binned.coords["sample_position"]
    upstream = binned.coords["source_position"]
    binned.coords["position"] = beam_frame(binned.coords["position"], sample_position, upstream)
    binned.coords["source_position"] = beam_frame(upstream, sample_position, upstream)
    binned.coords["sample_position"] = sc.vector([0.0, 0.0, 0.0], unit="m")
    return binned


def print_subpulse_table(pulses: sc.DataArray, events: sc.DataArray) -> None:
    """Print the incident energy, ki and time at the sample of each sub-pulse.

    Parameters
    ----------
    pulses:
        The sub-pulses of :func:`subpulse_timing`.
    events:
        The neutrons of :func:`.events_to_q_and_energy`.
    """
    print("  sub-pulse  t_sample [ms]  λi [Å]  Ei [meV]  ki [1/Å]  neutrons  intensity [n/s]")
    counts = np.bincount(events.coords["subpulse"].values, minlength=pulses.sizes["subpulse"])
    weights = np.bincount(
        events.coords["subpulse"].values, events.values, minlength=pulses.sizes["subpulse"]
    )
    for i, (t, wavelength, ei, n, w) in enumerate(
        zip(
            pulses.coords["time_at_sample"].to(unit="ms").values,
            pulses.coords["incident_wavelength"].to(unit="Å").values,
            pulses.coords["incident_energy"].to(unit="meV").values,
            counts,
            weights,
        )
    ):
        k = 2 * np.pi / wavelength
        print(f"  {i:9d}  {t:13.3f}  {wavelength:6.4f}  {ei:8.3f}  {k:8.4f}  {n:8d}  {w:15.4g}")


def _draw_plane(
    ax: Axes, events: sc.DataArray, edges: dict[str, sc.Variable], norm: LogNorm
) -> QuadMesh:
    """Draw the intensity of the neutrons in (Qx, Qz), summed over Qy, on ``edges``."""
    plane = events.hist(Qz=edges["Qz"], Qx=edges["Qx"])  # summed over Qy
    mesh = ax.pcolormesh(
        edges["Qx"].values,
        edges["Qz"].values,
        np.ma.masked_less_equal(plane.values, 0),  # empty bins in white
        cmap=plt.get_cmap("viridis").with_extremes(bad="white"),
        norm=norm,
    )
    ax.set(xlabel="Qx [1/Å]", ylabel="Qz [1/Å], along the beam", aspect="equal")
    ax.grid(True, alpha=0.3)
    return mesh


def plot_plane(
    events: sc.DataArray,
    pulses: sc.DataArray,
    title: str,
    step: float = Q_STEP,
    qy_half_width: float = QY_HALF_WIDTH,
) -> Figure:
    """Draw the intensity in (Qx, Qz) of every sub-pulse, summed over Qy and in a thin slice.

    The left panel sums over all Qy. The detector spans 25° below the beam to 15° above, and a
    neutron scattered out of the plane has less of kf in it, so each elastic circle spreads
    inward into a band. The right panel keeps ``|Qy| < qy_half_width`` only, where the neutrons
    of each sub-pulse lie on its own elastic circle, of radius ki through the origin, drawn in
    red.

    Both panels share one color scale, logarithmic, ``LOG_DECADES`` deep below the highest bin,
    in viridis with the empty bins in white, and their axes have the same scale, so that circles
    look round.

    Parameters
    ----------
    events:
        The neutrons of :func:`.events_to_q_and_energy`.
    pulses:
        The sub-pulses of :func:`subpulse_timing`, for the elastic circles.
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
    edges = {}
    for name in ("Qx", "Qz"):
        values = events.coords[name].to(unit="1/Å").values
        start, stop = np.floor(values.min() * 10) / 10, np.ceil(values.max() * 10) / 10
        edges[name] = bin_edges(name, start, stop, step, unit="1/Å")
    vmax = events.hist(Qz=edges["Qz"], Qx=edges["Qx"]).max().value
    norm = LogNorm(vmin=vmax / 10**LOG_DECADES, vmax=vmax, clip=True)
    qy = events.coords["Qy"].to(unit="1/Å").values
    (dim,) = events.dims
    in_slice = events[dim, np.flatnonzero(np.abs(qy) < qy_half_width)]

    fig, (summed_ax, slice_ax) = plt.subplots(
        1, 2, figsize=(16, 9), sharex=True, sharey=True, layout="constrained"
    )
    fig.suptitle(title)
    mesh = _draw_plane(summed_ax, events, edges, norm)
    _draw_plane(slice_ax, in_slice, edges, norm)
    summed_ax.set_title("Summed over Qy")
    slice_ax.set_title(f"|Qy| < {qy_half_width:g} 1/Å, with the elastic circle of each sub-pulse")
    angle = np.linspace(0.0, 2 * np.pi, 721)
    for wavelength in pulses.coords["incident_wavelength"].to(unit="Å").values:
        ki = 2 * np.pi / wavelength
        # Q = ki - kf with |kf| = ki: a circle of radius ki, centered on ki along the beam
        slice_ax.plot(ki * np.sin(angle), ki * (1 - np.cos(angle)), color="red", lw=0.4)
    slice_ax.set(xlim=summed_ax.get_xlim(), ylim=summed_ax.get_ylim())
    fig.colorbar(
        mesh, ax=[summed_ax, slice_ax], shrink=0.8, label=f"Intensity [n/s per ({step:g} 1/Å)²]"
    )
    return fig


def main(out: Path | None = None, ncount: float = 1e9, seed: int = 1) -> Path:
    """Run :class:`TRexWithLa2Ni7` in McStas, reduce its neutrons and save the (Qx, Qz) plane.

    The figure is ``trex_elastic_plane.png``, from :func:`plot_plane`.

    Parameters
    ----------
    out:
        Folder of the figure. If ``None``, the folder of the McStas run, next to its
        ``mccode.h5``.
    ncount:
        Number of McStas neutrons.
    seed:
        Seed of the random numbers.

    Returns
    -------
    :
        The path of the figure.
    """
    instrument = TRexWithLa2Ni7()  # the default control parameters
    name = "trex_elastic_plane"

    folder = Simulation(instrument).run(name, ncount=ncount, seed=seed).parent
    mcstas_data = ms.load_data(str(folder))
    bins = {monitor: monitor_bins(getattr(instrument, monitor)) for monitor in MONITORS}
    pulses = subpulse_timing(mcstas_data, bins, instrument.positions())

    events = events_to_q_and_energy(read_events_in_beam_frame(folder), pulses)
    print_subpulse_table(pulses, events)

    fig = plot_plane(
        events,
        pulses,
        f"T-Rex, La2Ni7: (Qx, Qz) of {pulses.sizes['subpulse']} sub-pulses\n"
        f"{instrument.running_conditions()}",
    )
    figure = output_folder(out, folder) / f"{name}.png"
    fig.savefig(figure, dpi=150)
    plt.close(fig)
    print(f"Saved {figure}")
    return figure


if __name__ == "__main__":
    main(out=parse_out(__doc__))
