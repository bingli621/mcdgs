"""The horizontal scattering plane of T-Rex with a La2Ni7 crystal, in (Qx, Qz).

:class:`TRexWithLa2Ni7` is :class:`.TRex` with the La2Ni7 crystal of :func:`.la2ni7` on
``sample_pos``, with the beamstop and the banana detector, which lists every neutron reaching
it. Its control parameters are those of :class:`.TRex`. Few neutrons reach a 5 mm crystal, so
McStas repeats each one that does ``SAMPLE_SPLIT`` times, with its weight divided among them.

The reduction follows :mod:`mcdgs.reduction`:

1. :func:`subpulse_timing` finds the sub-pulses at the pre-sample and beamstop monitors of McStas
   with :func:`.find_patches`, as :mod:`.trex_08_determine_ei` does, and gives each its incident
   energy, from its flight between the monitors with :func:`.incident_flight`, and the time it
   reaches the sample, with :func:`.time_at_sample`.
2. :func:`read_events` reads the neutrons of the banana detector, by pixel, in the beam frame of
   the sample: z along the beam between the pre-sample monitor and the sample, y up.
3. :func:`reduce_events` gives each pixel and sub-pulse the arrival window of
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
import mcstastox
import numpy as np
import scipp as sc
from matplotlib.axes import Axes
from matplotlib.collections import QuadMesh
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure
from mcstasscript.data.data import McStasData

from mcdgs.core.components import Single_crystal
from mcdgs.diagnostics.cli import output_folder, parse_out
from mcdgs.diagnostics.trex_07_neutrons_on_sample import find_patches
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.mcstas.simulation import Simulation
from mcdgs.reduction.binning import bin_edges
from mcdgs.reduction.incident_energy import incident_flight
from mcdgs.reduction.inelastic import (
    DEFAULT_ENERGY_TRANSFER_RATIO,
    QE_GRAPH,
    arrival_window,
    assign_subpulse,
    time_at_sample,
)
from mcdgs.samples import la2ni7
from mcdgs.tof.reader import monitor_bins

MONITORS = ("pre_sample_monitor", "beamstop_monitor")  # in beam order
DETECTOR = "banana_detector"
SAMPLE = "sample_pos"
MIN_FRACTION = 0.02  # a patch weaker than this fraction of the strongest is not a sub-pulse
Q_STEP = 0.02  # 1/Å, the bin width in Qx and Qz
QY_HALF_WIDTH = 0.05  # 1/Å, half the width of the slice in Qy near the horizontal plane
LOG_DECADES = 4  # the log color scale spans this many decades below the highest bin
SAMPLE_SPLIT = 100  # McStas SPLIT on the crystal: each neutron reaching it is repeated this often


class TRexWithLa2Ni7(TRex):
    """T-Rex with a La2Ni7 crystal, from the source to the banana detector.

    Parameters
    ----------
    sample:
        The crystal. If ``None``, :func:`.la2ni7`, split ``SAMPLE_SPLIT`` times.
    params:
        Control parameters.
    """

    name = "trex_la2ni7"  # the McStas instrument name

    def __init__(
        self, sample: Single_crystal | None = None, params: Params = DEFAULT_PARAMS
    ) -> None:
        if sample is None:
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
    the intensity of the strongest; the n-th at one monitor is the n-th at the other, as in
    :func:`.trex_08_determine_ei.subpulse_flights`.

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
    found = []
    for name in MONITORS:
        # McStasScript has the map as (wavelength, time), the y axis first
        intensity = np.asarray(ms.name_search(name, mcstas_data).Intensity, dtype=float).T
        patches = find_patches(intensity, bins[name])
        found.append((patches, patches.intensity >= min_fraction * patches.intensity.max()))
    (first, keep_first), (second, keep_second) = found
    if keep_first.sum() != keep_second.sum():
        raise ValueError(
            f"{keep_first.sum()} sub-pulses at the pre-sample monitor, "
            f"{keep_second.sum()} at the beamstop monitor"
        )
    time_at_pre_sample = sc.array(dims=["subpulse"], values=first.time[keep_first], unit="ms")
    pre_sample_position = sc.vector(positions["pre_sample_monitor"], unit="m")
    flight = incident_flight(
        pre_sample_position,
        time_at_pre_sample,
        sc.vector(positions["beamstop_monitor"], unit="m"),
        sc.array(dims=["subpulse"], values=second.time[keep_second], unit="ms"),
    )
    if (flight["tof"] <= sc.scalar(0.0, unit="ms")).any().value:
        raise ValueError("A sub-pulse reaches the beamstop monitor first")
    return sc.DataArray(
        sc.array(dims=["subpulse"], values=first.intensity[keep_first], unit="1/s"),  # n/s
        coords={
            **flight,
            "time_at_pre_sample_monitor": time_at_pre_sample,
            "time_at_sample": time_at_sample(
                time_at_pre_sample,
                pre_sample_position,
                sc.vector(positions[SAMPLE], unit="m"),
                flight["incident_energy"],
            ),
        },
    )


def beam_frame(
    position: sc.Variable, sample_position: sc.Variable, upstream_position: sc.Variable
) -> sc.Variable:
    """Return positions in the beam frame of the sample.

    Its origin is the sample, z points along the beam, from ``upstream_position`` to the sample,
    y up, as the y of the lab, and x to the left of the beam, ``y × z``, as in McStas.

    Parameters
    ----------
    position:
        The positions, as vectors in the lab frame.
    sample_position:
        Where the sample is.
    upstream_position:
        A point on the beam before the sample, such as a monitor.

    Returns
    -------
    :
        The positions in the beam frame, in the unit of ``position``.

    Examples
    --------
    >>> local = beam_frame(
    ...     sc.vector([2.0, 0.0, 1.0], unit="m"),
    ...     sc.vector([1.0, 0.0, 1.0], unit="m"),
    ...     sc.vector([0.0, 0.0, 0.0], unit="m"),
    ... )
    >>> np.round(local.value, 4).tolist()  # 45° to the left of a beam at 45°
    [0.7071, 0.0, 0.7071]
    """
    beam = sample_position - upstream_position.to(unit=sample_position.unit)
    z = beam / sc.norm(beam)
    y = sc.vector([0.0, 1.0, 0.0])
    x = sc.cross(y, z)
    x = x / sc.norm(x)
    y = sc.cross(z, x)  # exactly across the beam, should it not be horizontal
    relative = position - sample_position.to(unit=position.unit)
    return sc.spatial.as_vectors(sc.dot(relative, x), sc.dot(relative, y), sc.dot(relative, z))


def read_events(folder: str | Path) -> sc.DataArray:
    """Read the neutrons of the banana detector of a McStas run, in the beam frame.

    Parameters
    ----------
    folder:
        Run folder with the ``mccode.h5`` to read.

    Returns
    -------
    :
        The neutrons, binned by ``pixel_id``, each with its weight in n/s and its time
        ``time_at_detector``. The pixels have their ``position``, and the data its
        ``sample_position`` and ``source_position``: the pre-sample monitor, on the beam. Times in
        ms, positions in
        the frame of :func:`beam_frame`.
    """
    with mcstastox.Read(str(folder)) as data:
        # the pre-sample monitor as the source: the beam from it to the sample is the incident
        # beam, which the curved guide does not give from the source itself
        binned = data.export_scipp(MONITORS[0], SAMPLE, component_name=DETECTOR)["events"]
    sample_position = binned.coords["sample_position"]
    upstream = binned.coords["source_position"]
    binned.coords["position"] = beam_frame(binned.coords["position"], sample_position, upstream)
    binned.coords["source_position"] = beam_frame(upstream, sample_position, upstream)
    binned.coords["sample_position"] = sc.vector([0.0, 0.0, 0.0], unit="m")
    # in ms, as the times of the sub-pulses at the sample
    binned.bins.coords["time_at_detector"] = binned.bins.coords.pop("t").to(unit="ms")
    binned.bins.unit = "1/s"  # McStas weights are n/s
    return binned


def reduce_events(
    binned: sc.DataArray,
    pulses: sc.DataArray,
    energy_transfer_ratio: tuple[float, float] = DEFAULT_ENERGY_TRANSFER_RATIO,
) -> sc.DataArray:
    """Keep the neutrons inside the arrival window of a sub-pulse, with their Q and E.

    Parameters
    ----------
    binned:
        The neutrons of :func:`read_events`.
    pulses:
        The sub-pulses of :func:`subpulse_timing`.
    energy_transfer_ratio:
        The smallest and the largest ``E / Ei`` of the arrival windows, as
        :func:`.arrival_window`.

    Returns
    -------
    :
        The neutrons inside a window, along ``event``, with ``subpulse``, ``incident_energy``,
        ``energy_transfer``, ``Q``, ``Qx``, ``Qy`` and ``Qz``, in the beam frame.
    """
    L2 = sc.norm(binned.coords["position"] - binned.coords["sample_position"])
    window = arrival_window(
        pulses.coords["time_at_sample"],
        L2,
        pulses.coords["incident_energy"],
        energy_transfer_ratio=energy_transfer_ratio,
    )
    # every neutron gets the id and the position of its pixel, then they are one list
    binned = binned.copy()
    for name in ("pixel_id", "position"):
        binned.bins.coords[name] = sc.bins_like(binned, binned.coords[name])
    events = binned.bins.concat().value
    (dim,) = events.dims
    assigned = assign_subpulse(
        events.coords["pixel_id"],
        events.coords["time_at_detector"],
        binned.coords["pixel_id"],
        window["earliest_arrival"],
        window["latest_arrival"],
    )
    inside = assigned["inside"].values
    events = events[dim, np.flatnonzero(inside)].copy()
    subpulse = assigned["subpulse"].values[inside]
    events.coords["sample_position"] = binned.coords["sample_position"]
    events.coords["source_position"] = binned.coords["source_position"]
    events.coords["subpulse"] = sc.array(dims=[dim], values=subpulse)
    for name in ("incident_energy", "time_at_sample"):
        coord = pulses.coords[name]
        events.coords[name] = sc.array(dims=[dim], values=coord.values[subpulse], unit=coord.unit)
    return events.transform_coords(
        ["energy_transfer", "Q", "Qx", "Qy", "Qz"], graph=QE_GRAPH, keep_intermediate=False
    )


def print_subpulse_table(pulses: sc.DataArray, events: sc.DataArray) -> None:
    """Print the incident energy, ki and time at the sample of each sub-pulse.

    Parameters
    ----------
    pulses:
        The sub-pulses of :func:`subpulse_timing`.
    events:
        The neutrons of :func:`reduce_events`.
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
        The neutrons of :func:`reduce_events`.
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

    events = reduce_events(read_events(folder), pulses)
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
