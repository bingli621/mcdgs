"""The neutrons that reach the sample of T-Rex, from McStas and from tof.

:class:`TRexToSample` is :class:`.TRex` up to its sample monitor, at ``sample_pos``: without a
sample, and without the beamstop and the banana detector after it. Its control parameters are
those of :class:`.TRex`.

tof has no guides: every neutron its choppers let through reaches the sample. As in
:mod:`.trex_04_beam_monitor`, its counts are scaled at the source monitor, to the total of
McStas there, where both see the same beam. At the sample monitor, a transmission on top of that
scale stands for what the guides take out of the beam. It is tuned so that tof matches McStas:
by default, so that both have the same total. The transmission is then that of the guides in
McStas, and the shape in time and wavelength that of each model.

The choppers cut the beam into sub-pulses, a bin or two wide and about 4 ms apart.
:func:`.plot_comparison` draws both models twice, with the points of the wavelength spectrum and
the time profile joined by lines: first all on linear scales, then all on log scales, the maps
``LOG_DECADES`` deep, so that the weak bins show too. :func:`.find_patches` finds the patches of
each model, the sub-pulses and the weak ones between them, and :func:`plot_patches` draws the
center of mass of each in time and wavelength, with its width as error bars.

:func:`subpulse_windows` reduces the tof run as :mod:`.trex_09_elastic_plane` reduces McStas:
the sub-pulses at the pre-sample and beamstop monitors give their incident energy and the time
they reach the sample, and :func:`.arrival_window` when the neutrons of each reach the banana
detector. :func:`windows_table` lists them.

:func:`main` runs both models, reads the intensity of each on the bins of the monitors, scales
tof at the source monitor, tunes the transmission at the sample monitor and saves the figures.
Run the module to do it::

    pixi run python -m mcdgs.diagnostics.trex_07_neutrons_on_sample
"""

from pathlib import Path

import matplotlib.pyplot as plt
import mcstasscript as ms
import scipp as sc
import tof
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from mcdgs.diagnostics.cli import output_folder, parse_out
from mcdgs.diagnostics.compare import MonitorComparison, plot_comparison, read_scaled_comparisons
from mcdgs.instruments.trex.detectors import BANANA_RADIUS
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.mcstas.simulation import Simulation
from mcdgs.reduction.inelastic import DEFAULT_ENERGY_TRANSFER_RATIO, arrival_window, time_at_sample
from mcdgs.reduction.patches import find_patches, subpulse_flight
from mcdgs.tof.builder import build_tof
from mcdgs.tof.reader import monitor_bins, read_tof_monitor

LOG_DECADES = 3  # the log color scale of the maps spans this many decades below the highest bin
MONITORS = ("pre_sample_monitor", "beamstop_monitor")  # give the incident energy, in beam order


# --8<-- [start:trex_to_sample]
class TRexToSample(TRex):
    """T-Rex from the source to the sample monitor, without a sample.

    The beamstop and the banana detector, after the sample, are left out.

    Parameters
    ----------
    params:
        Control parameters.
    """

    name = "trex_to_sample"  # the McStas instrument name

    def __init__(self, params: Params = DEFAULT_PARAMS) -> None:
        super().__init__(params=params)
        # after the sample monitor
        del self.beamstop
        del self.banana_detector
        # --8<-- [end:trex_to_sample]


def plot_patches(comparison: MonitorComparison, title: str) -> Figure:
    """Draw the center of mass of each patch of both models, in time and wavelength.

    Each patch is a point at its mean time and wavelength, with its width in each as gray error
    bars, from :func:`.find_patches`. Its color is its intensity, on one log scale for both models:
    McStas has small filled circles, tof larger hollow squares with the color as their outline, so
    that where both have a patch the square frames the circle.

    Parameters
    ----------
    comparison:
        Both maps.
    title:
        The title of the figure.

    Returns
    -------
    :
        The figure.
    """
    mcstas = find_patches(comparison.mcstas_intensity, comparison.bins)
    tof = find_patches(comparison.tof_intensity, comparison.bins)
    norm = LogNorm(
        vmin=min(mcstas.intensity.min(), tof.intensity.min()),
        vmax=max(mcstas.intensity.max(), tof.intensity.max()),
    )
    cmap = plt.get_cmap("turbo")

    fig, ax = plt.subplots(figsize=(11, 6), layout="constrained")
    for p in (mcstas, tof):
        ax.errorbar(
            p.time,
            p.wavelength,
            xerr=p.time_width,
            yerr=p.wavelength_width,
            fmt="none",
            ecolor="0.6",
            capsize=2,
            zorder=1,
        )
    # McStas: small filled circles; tof: larger hollow squares, the color as the outline
    points = ax.scatter(
        mcstas.time,
        mcstas.wavelength,
        c=mcstas.intensity,
        cmap=cmap,
        norm=norm,
        marker="o",
        s=25,
        edgecolors="black",
        linewidths=0.3,
        zorder=2,
    )
    ax.scatter(
        tof.time,
        tof.wavelength,
        marker="s",
        s=110,
        facecolors="none",
        edgecolors=cmap(norm(tof.intensity)),
        linewidths=1.8,
        zorder=3,
    )
    # the legend shows the markers only, not the color of one patch
    handles = [
        Line2D(
            [],
            [],
            linestyle="none",
            marker="o",
            markersize=5,
            markerfacecolor="white",
            markeredgecolor="black",
            label=f"McStas, {mcstas.time.size} patches",
        ),
        Line2D(
            [],
            [],
            linestyle="none",
            marker="s",
            markersize=10,
            markerfacecolor="none",
            markeredgecolor="black",
            markeredgewidth=1.8,
            label=f"tof, {tof.time.size} patches",
        ),
    ]
    fig.colorbar(points, ax=ax, label="Intensity of the patch [n/s]")
    ax.set(title=title, xlabel="Time [ms]", ylabel="Wavelength [Å]")
    ax.grid(True, alpha=0.4)
    ax.legend(handles=handles)
    return fig


def subpulse_windows(
    instrument: TRex,
    tof_result: tof.Result,
    energy_transfer_ratio: tuple[float, float] = DEFAULT_ENERGY_TRANSFER_RATIO,
) -> sc.DataArray:
    """Give each sub-pulse of a tof run its arrival window at the banana detector.

    As :mod:`.trex_09_elastic_plane` does with McStas: :func:`.find_patches` finds the
    sub-pulses at both monitors, in the first pulse, :func:`.subpulse_flight` gives their
    incident energy from the flight between them, :func:`.time_at_sample` the time they reach
    the sample, and :func:`.arrival_window` when their neutrons reach the detector,
    ``BANANA_RADIUS`` from the sample.

    Parameters
    ----------
    instrument:
        The instrument of the run, with the monitors of ``MONITORS``.
    tof_result:
        The tof run.
    energy_transfer_ratio:
        The smallest and the largest ``E / Ei`` of the windows, as :func:`.arrival_window`.

    Returns
    -------
    :
        The sub-pulses along ``subpulse``, with the coordinates of :func:`.subpulse_flight`,
        ``time_at_sample``, and ``earliest_arrival`` and ``latest_arrival`` at the detector, in
        ms. The data are tof neutrons at the pre-sample monitor, not n/s.
    """
    positions = {name: sc.vector(xyz, unit="m") for name, xyz in instrument.positions().items()}
    start, end = (
        find_patches(read_tof_monitor(tof_result[name], bins), bins)
        for name in MONITORS
        for bins in [monitor_bins(getattr(instrument, name))]
    )
    pulses = subpulse_flight(start, positions[MONITORS[0]], end, positions[MONITORS[1]])
    pulses.unit = None  # tof neutrons, not n/s
    pulses.coords["time_at_sample"] = time_at_sample(
        pulses.coords["time_at_start"],
        positions[MONITORS[0]],
        positions["sample_pos"],
        pulses.coords["incident_energy"],
    )
    window = arrival_window(
        pulses.coords["time_at_sample"],
        sc.scalar(BANANA_RADIUS, unit="m"),
        pulses.coords["incident_energy"],
        energy_transfer_ratio=energy_transfer_ratio,
    )
    pulses.coords.update(window)
    return pulses


def windows_table(windows: sc.DataArray) -> str:
    """Return the arrival windows of :func:`subpulse_windows` as a Markdown table.

    Each row is a sub-pulse: its incident wavelength and energy, when it reaches the sample,
    and when its window at the detector opens and closes. The last column is the time from
    the close of its window to the opening of the next: negative where the two overlap.

    Parameters
    ----------
    windows:
        The sub-pulses of :func:`subpulse_windows`.

    Returns
    -------
    :
        The table.
    """
    wavelength = windows.coords["incident_wavelength"].to(unit="Å").values
    energy = windows.coords["incident_energy"].to(unit="meV").values
    at_sample = windows.coords["time_at_sample"].to(unit="ms").values
    opens = windows.coords["earliest_arrival"].to(unit="ms").values
    closes = windows.coords["latest_arrival"].to(unit="ms").values
    rows = [
        "| Sub-pulse | λi | Ei | At the sample | Window opens | Window closes | To the next |",
        "| --------- | -- | -- | ------------- | ------------ | ------------- | ----------- |",
    ]
    for i in range(windows.sizes["subpulse"]):
        gap = ""
        if i + 1 < len(opens):
            # a minus sign, as in the docs
            gap = f"{opens[i + 1] - closes[i]:.3f} ms".replace("-", "\u2212")
        rows.append(
            f"| {i} | {wavelength[i]:.4f} Å | {energy[i]:.3f} meV | {at_sample[i]:.3f} ms "
            f"| {opens[i]:.3f} ms | {closes[i]:.3f} ms | {gap} |"
        )
    return "\n".join(rows) + "\n"


def main(
    out: Path | None = None,
    ncount: float = 1e8,
    neutrons: int = 10_000_000,
    seed: int = 1,
    transmission: float | None = None,
) -> list[Path]:
    """Run both models of :class:`TRexToSample`, print the transmission and save the figures.

    The figures are ``trex_neutrons_on_sample.png``, on linear scales, and
    ``trex_neutrons_on_sample_log.png``, on log scales, from :func:`.plot_comparison`; and
    ``trex_neutrons_on_sample_patches.png``, the center of mass of each patch, from
    :func:`plot_patches`. ``trex_neutrons_on_sample_windows.txt`` is the table of the arrival
    windows of tof, from :func:`windows_table`.

    Parameters
    ----------
    out:
        Folder of the figures. If ``None``, the folder of the McStas run, next to its
        ``mccode.h5``.
    ncount:
        Number of McStas neutrons. Many: the guides and the choppers let through a small
        fraction of the beam.
    neutrons:
        Number of tof neutrons.
    seed:
        Seed of the random numbers of both models.
    transmission:
        The factor on tof at the sample monitor, on top of the scale of the source monitor: the
        fraction of the beam the guides let through, which tof does not model. If ``None``,
        tuned so that tof has the total of McStas at the sample monitor.

    Returns
    -------
    :
        The paths of the figures: linear, log, then the patches; and of the table of the
        arrival windows.
    """
    instrument = TRexToSample()  # the default control parameters
    name = "trex_neutrons_on_sample"

    # run McStas and load its monitors, in n/s per bin
    folder = Simulation(instrument).run(name, ncount=ncount, seed=seed).parent
    mcstas_data = ms.load_data(str(folder))

    # run tof, then read the monitors of both on the bins of the McStas maps: time, wavelength
    tof_result = build_tof(instrument, neutrons=neutrons, seed=seed).run()
    # tof's counts are scaled to the total of McStas at the source monitor, in n/s; at the
    # sample monitor, the same scale, and the transmission on top of it
    comparisons = read_scaled_comparisons(
        mcstas_data,
        tof_result,
        {
            "source_monitor": monitor_bins(instrument.source_monitor),
            "sample_monitor": monitor_bins(instrument.sample_monitor),
        },
        ncount,
        neutrons,
        transmission=transmission,
    )
    at_source, comparison = comparisons["source_monitor"], comparisons["sample_monitor"]
    transmission = comparison.transmission

    # print the totals, and the scale and transmission of tof
    mcstas_total = comparison.mcstas_intensity.sum()
    tof_counts = comparison.tof_intensity.sum() / comparison.scale / transmission
    print(
        f"McStas total {mcstas_total:.4g} n/s, "
        f"{mcstas_total / at_source.mcstas_intensity.sum():.4g} of that at the source; "
        f"tof {tof_counts:.0f} of {neutrons:,} neutrons, "
        f"{comparison.scale:.4g} n/s each, transmission {transmission:.4g}"
    )

    title = f"T-Rex flux at the sample monitor\n{instrument.running_conditions()}"
    figures = []
    for suffix, log_decades, log_y in (("", None, False), ("_log", LOG_DECADES, True)):
        figure = output_folder(out, folder) / f"{name}{suffix}.png"
        fig = plot_comparison(comparison, title, log_decades=log_decades, connect=True, log_y=log_y)
        fig.savefig(figure, dpi=150)
        plt.close(fig)
        print(f"Saved {figure}")
        figures.append(figure)

    figure = output_folder(out, folder) / f"{name}_patches.png"
    fig = plot_patches(
        comparison, f"Center of mass of each patch\n{instrument.running_conditions()}"
    )
    fig.savefig(figure, dpi=150)
    plt.close(fig)
    print(f"Saved {figure}")
    figures.append(figure)

    table = output_folder(out, folder) / f"{name}_windows.txt"
    table.write_text(windows_table(subpulse_windows(instrument, tof_result)))
    print(table.read_text(), end="")
    print(f"Saved {table}")
    figures.append(table)
    return figures


if __name__ == "__main__":
    main(out=parse_out(__doc__))
