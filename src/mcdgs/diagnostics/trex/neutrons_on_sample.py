"""The neutrons that reach the sample of T-Rex, from McStas and from tof.

:class:`TRexToSample` is :class:`.TRex` up to its sample monitor, at ``sample_pos``: without a
sample, and without the beamstop and the banana detector after it. Its control parameters are
those of :class:`.TRex`. After the sample monitor it has :func:`.sample_xy_monitor`, in the same
place, which records where the neutrons cross the sample position: one map of x against y.
:func:`.read_xy_map` reads that map from the McStas run and :func:`.plot_xy_map` draws it, with
its profile along each axis. tof has no x and y, so the map is that of McStas alone.

tof has no guides: every neutron its choppers let through reaches the sample. As in
:mod:`.diagnostics.trex.beam_monitor`, its counts are scaled at the source monitor, to the total of
McStas there, where both see the same beam. At the sample monitor, a transmission on top of that
scale stands for what the guides take out of the beam. It is tuned so that tof matches McStas:
by default, so that both have the same total. The transmission is then that of the guides in
McStas, and the shape in time and wavelength that of each model.

The choppers cut the beam into sub-pulses, a bin or two wide and about 4 ms apart.
:func:`.plot_comparison` draws both models twice, with the points of the wavelength spectrum and
the time profile joined by lines: first all on linear scales, then all on log scales, the maps
``LOG_DECADES`` deep, so that the weak bins show too. :func:`.find_patches` finds the patches of
each model, the sub-pulses and any weak ones between them, and :func:`.plot_patches` draws the
center of mass of each in time and wavelength, with its width as error bars.

:func:`subpulse_windows` reduces the tof run as :mod:`.diagnostics.trex.diffraction` reduces McStas:
the sub-pulses at the pre-sample and beamstop monitors give their incident energy and the time
they reach the sample, and :func:`.arrival_window` when the neutrons of each reach the banana
detector. :func:`.windows_table` lists them.

:func:`main` runs both models, reads the intensity of each on the bins of the monitors, scales
tof at the source monitor, tunes the transmission at the sample monitor and saves the figures,
that of the x and y map among them.
Run the module to do it::

    pixi run python -m mcdgs.diagnostics.trex.neutrons_on_sample
"""

from pathlib import Path

import scipp as sc
import tof

from mcdgs.diagnostics import cli, compare
from mcdgs.diagnostics.neutrons_on_sample import (
    elastic_arrival,
    plot_patches,
    plot_xy_map,
    read_xy_map,
    sample_xy_monitor,
    windows_table,
)
from mcdgs.instruments.trex.detectors import BANANA_RADIUS
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.instruments.trex.monitors import SAMPLE_MONITOR_HEIGHT, SAMPLE_MONITOR_WIDTH
from mcdgs.mcstas import Simulation, read_mcstas_maps
from mcdgs.reduction.binning import monitor_bins
from mcdgs.reduction.inelastic import DEFAULT_ENERGY_TRANSFER_RATIO, arrival_window, time_at_sample
from mcdgs.reduction.patches import find_patches, subpulse_flight
from mcdgs.tof import build_tof, read_tof_maps

LOG_DECADES = 3  # the log color scale of the maps spans this many decades below the highest bin
MONITORS = ("pre_sample_monitor", "beamstop_monitor")  # give the incident energy, in beam order


# --8<-- [start:instrument]
class TRexToSample(TRex):
    """T-Rex from the source to the sample monitor, without a sample.

    The beamstop and the banana detector, after the sample, are left out. In their place is
    :func:`.sample_xy_monitor`, at ``sample_pos`` as the sample monitor.

    Parameters
    ----------
    params:
        Control parameters.
    source_masked:
        As in :class:`.TRex`.
    """

    name = "trex_to_sample"  # the McStas instrument name

    def __init__(self, params: Params = DEFAULT_PARAMS, source_masked: bool = False) -> None:
        super().__init__(params=params, source_masked=source_masked)
        # after the sample monitor
        del self.beamstop
        del self.banana_detector
        self.sample_xy_monitor = sample_xy_monitor(SAMPLE_MONITOR_WIDTH, SAMPLE_MONITOR_HEIGHT)


def define_instrument() -> TRexToSample:
    """Return the instrument of this diagnostic: :class:`TRexToSample` at its defaults."""
    params = Params()  # the default control parameters
    return TRexToSample(params=params)


# --8<-- [end:instrument]


def subpulse_windows(
    instrument: TRex,
    tof_run: tof.Result,
    energy_transfer_ratio: tuple[float, float] = DEFAULT_ENERGY_TRANSFER_RATIO,
) -> sc.DataArray:
    """Give each sub-pulse of a tof run its arrival window at the banana detector.

    As :mod:`.diagnostics.trex.diffraction` does with McStas: :func:`.find_patches` finds the
    sub-pulses at both monitors, in the first pulse, :func:`.subpulse_flight` gives their
    incident energy from the flight between them, :func:`.time_at_sample` the time they reach
    the sample, and :func:`.arrival_window` when their neutrons reach the detector,
    ``BANANA_RADIUS`` from the sample.

    Parameters
    ----------
    instrument:
        The instrument of the run, with the monitors of ``MONITORS``.
    tof_run:
        The tof run.
    energy_transfer_ratio:
        The smallest and the largest ``E / Ei`` of the windows, as :func:`.arrival_window`.

    Returns
    -------
    :
        The sub-pulses along ``subpulse``, with the coordinates of :func:`.subpulse_flight`,
        ``time_at_sample``, and ``earliest_arrival``, ``elastic_arrival``, of a neutron that
        keeps its speed, and ``latest_arrival`` at the detector, in ms. The data are tof
        neutrons at the pre-sample monitor, not n/s.
    """
    positions = {name: sc.vector(xyz, unit="m") for name, xyz in instrument.positions().items()}
    bins = {name: monitor_bins(getattr(instrument, name).bin_limits) for name in MONITORS}
    tof_maps = read_tof_maps(tof_run, bins)
    start, end = (find_patches(tof_maps[name].values, bins[name]) for name in MONITORS)
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
    # a neutron that flies through the sample without changing its speed
    pulses.coords["elastic_arrival"] = elastic_arrival(
        pulses.coords["time_at_sample"],
        pulses.coords["incident_wavelength"],
        sc.scalar(BANANA_RADIUS, unit="m"),
    )
    return pulses


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
    :func:`.plot_patches`. ``trex_neutrons_on_sample_windows.txt`` is the table of the arrival
    windows of tof, from :func:`.windows_table`. ``trex_neutrons_on_sample_xy.png`` is the x and y
    map of McStas at the sample position, from :func:`.plot_xy_map`.

    Parameters
    ----------
    out:
        Folder of the figures. If ``None``, ``trex_neutrons_on_sample`` in the results folder of
        :meth:`.Workspace.find`.
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
        The paths of the figures: linear, log, then the patches; of the table of the
        arrival windows; and of the x and y map.
    """
    instrument = define_instrument()
    name = "trex_neutrons_on_sample"

    # run McStas
    sim = Simulation(instrument)
    mcstas_run = sim.run(folder_name=name, ncount=ncount, seed=seed)

    # run tof, then read the monitors of both on the bins of the McStas maps: time, wavelength
    tof_model = build_tof(instrument, neutrons=neutrons, seed=seed)
    tof_run = tof_model.run()
    # tof's counts are scaled to the total of McStas at the source monitor, in n/s; at the
    # sample monitor, the same scale, and the transmission on top of it
    bins = {
        name: monitor_bins(getattr(instrument, name).bin_limits)
        for name in ("source_monitor", "sample_monitor")
    }
    mcstas_maps = read_mcstas_maps(mcstas_run, bins)
    tof_maps = read_tof_maps(tof_run, bins)
    mcstas_source, tof_source = mcstas_maps["source_monitor"], tof_maps["source_monitor"]
    scale = compare.rescale_tof_to_mcstas(tof_source, mcstas_source)
    at_source = compare.mcstas_vs_tof(
        "source_monitor", mcstas_source, tof_source, ncount, neutrons, scale
    )
    mcstas_sample, tof_sample = mcstas_maps["sample_monitor"], tof_maps["sample_monitor"]
    if transmission is None:  # tuned: the transmission that gives tof the total of McStas
        transmission = compare.rescale_tof_to_mcstas(tof_sample, mcstas_sample) / scale
    comparison = compare.mcstas_vs_tof(
        "sample_monitor", mcstas_sample, tof_sample, ncount, neutrons, scale, transmission
    )

    # print the totals, and the scale and transmission of tof
    mcstas_total = comparison.first.values.sum()
    tof_counts = comparison.second.values.sum() / scale / transmission
    print(
        f"McStas total {mcstas_total:.4g} n/s, "
        f"{mcstas_total / at_source.first.values.sum():.4g} of that at the source; "
        f"tof {tof_counts:.0f} of {neutrons:,} neutrons, "
        f"{scale:.4g} n/s each, transmission {transmission:.4g}"
    )

    folder = cli.output_folder(out, name)
    conditions = instrument.running_conditions()
    files = [
        cli.save_figure(
            compare.plot_comparison(
                comparison,
                f"T-Rex flux at the sample monitor\n{conditions}",
                log_decades=log_decades,
                connect=True,
                log_y=log_y,
            ),
            folder / f"{name}{suffix}.png",
        )
        for suffix, log_decades, log_y in (("", None, False), ("_log", LOG_DECADES, True))
    ]

    windows = subpulse_windows(instrument, tof_run)
    fig = plot_patches(comparison, f"Center of mass of each patch\n{conditions}", windows=windows)
    files.append(cli.save_figure(fig, folder / f"{name}_patches.png"))

    table = windows_table(windows)
    print(table, end="")
    files.append(cli.save_text(table, folder / f"{name}_windows.txt"))

    xy_map = read_xy_map(instrument.sample_xy_monitor, mcstas_run)
    fig = plot_xy_map(xy_map, f"T-Rex beam at the sample position, McStas\n{conditions}")
    files.append(cli.save_figure(fig, folder / f"{name}_xy.png"))
    return files


if __name__ == "__main__":
    main(out=cli.parse_out(__doc__))
