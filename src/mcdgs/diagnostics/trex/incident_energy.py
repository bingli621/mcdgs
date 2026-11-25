"""The incident energy of T-Rex, from the direct beam at the pre-sample and beamstop monitors.

:class:`TRexToBeamstopMonitor` is :class:`.TRex` with the La2Ni7 crystal of :func:`.la2ni7` on
``sample_pos``, without the beamstop and the banana detector: the beamstop monitor is the last
monitor that counts. Its control parameters are those of :class:`.TRex`.

The beamstop monitor is ``BEAMSTOP_DISTANCE`` after the sample in space, but before it in beam
order: it counts the direct beam as if there were no sample. So both monitors see the same
sub-pulses, the beamstop monitor later by the flight over the distance between them.

As in :mod:`.diagnostics.trex.neutrons_on_sample`, tof's counts are scaled at the source monitor, to the
total of McStas there. A transmission on top of that scale stands for what the guides take out of
the beam, which tof does not model: by default, it is tuned so that tof has the total of McStas at
the pre-sample monitor, and the same factor holds at the beamstop monitor.

The incident energy of each sub-pulse follows from its flight between the monitors, and also
from its flight from the source to the pre-sample monitor alone, leaving the source at
``t0_offset``, the time the choppers are phased to.
:func:`subpulse_flights` finds the sub-pulses at each monitor with :func:`.find_patches` and pairs
them with :func:`.subpulse_flight`, both of :mod:`mcdgs.reduction.patches`.
:func:`.incident_flight`, of :mod:`mcdgs.reduction.incident_energy`, turns the start and end of
each flight, and the times there, into the flight path, the time of flight, the incident
wavelength and the incident energy, through a graph of scippneutron's kernels. The mean
wavelength of each sub-pulse at the pre-sample monitor, which the monitors record, gives the
incident energy it should have.

:func:`.plot_time_profiles` draws the intensity of both monitors against time,
:func:`.plot_incident_energies` the incident energy of each sub-pulse, and :func:`.save_graph` the
graph, all of :mod:`mcdgs.diagnostics.incident_energy`. Run the module to do it::

    pixi run python -m mcdgs.diagnostics.trex.incident_energy
"""

from pathlib import Path

import scipp as sc
from scippneutron.conversion import tof as tof_kernels

from mcdgs.diagnostics import cli, compare
from mcdgs.diagnostics.incident_energy import plot_incident_energies, plot_time_profiles, save_graph
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.mcstas import Simulation, read_mcstas_maps
from mcdgs.reduction.binning import monitor_bins
from mcdgs.reduction.incident_energy import incident_flight
from mcdgs.reduction.patches import MIN_FRACTION, find_patches, subpulse_flight
from mcdgs.samples import la2ni7
from mcdgs.tof import build_tof, read_tof_maps

MONITORS = ("pre_sample_monitor", "beamstop_monitor")  # in beam order
MONITOR_LABELS = {
    "pre_sample_monitor": "Pre-sample monitor",
    "beamstop_monitor": "Beamstop monitor",
}
# the methods of the incident energy, the prefix of the coordinates of each flight
TWO_MONITORS, SINGLE_MONITOR = "from_two_monitors", "from_single_monitor"
METHODS = {  # each method, its color and its label
    TWO_MONITORS: ("tab:blue", "between the monitors"),
    SINGLE_MONITOR: ("tab:orange", "pre-sample monitor alone"),
}
FLIGHT_COORDS = ("Ltotal", "tof", "incident_wavelength", "incident_energy")  # of each flight


# --8<-- [start:instrument]
class TRexToBeamstopMonitor(TRex):
    """T-Rex with a La2Ni7 crystal, from the source to the beamstop monitor.

    The beamstop and the banana detector, which only see what the sample scatters or lets
    through, are left out.

    The crystal is that of :func:`.la2ni7`.

    Parameters
    ----------
    params:
        Control parameters.
    source_masked:
        As in :class:`.TRex`. True by default: the masked source of
        :func:`.sources.masked_source` spends the rays where the choppers let neutrons through.
    """

    name = "trex_to_beamstop_monitor"  # the McStas instrument name

    def __init__(self, params: Params = DEFAULT_PARAMS, source_masked: bool = True) -> None:
        super().__init__(sample=la2ni7(name="sample"), params=params, source_masked=source_masked)
        # after the sample: the beamstop monitor is the last one that counts
        del self.beamstop
        del self.banana_detector


def define_instrument() -> TRexToBeamstopMonitor:
    """Return the instrument of this diagnostic: :class:`TRexToBeamstopMonitor` at its defaults."""
    params = Params()  # the default control parameters
    return TRexToBeamstopMonitor(params=params)


# --8<-- [end:instrument]


def subpulse_flights(
    pre_sample: compare.MonitorComparison,
    beamstop: compare.MonitorComparison,
    positions: dict[str, list[float]],
    time_on_source: sc.Variable,
    model: str,
    min_fraction: float = MIN_FRACTION,
) -> sc.DataArray:
    """Pair the sub-pulses of one model at both monitors, and compute their incident energy.

    The sub-pulses are the patches of :func:`.find_patches` with at least ``min_fraction`` of
    the intensity of the strongest, paired by :func:`.subpulse_flight`: the flight between the
    monitors is shorter than the time between sub-pulses, so the n-th sub-pulse at one is the
    n-th at the other.

    The incident energy comes from two flights, each through :func:`.incident_flight`: from the
    pre-sample to the beamstop monitor, and from the source to the pre-sample monitor alone,
    every sub-pulse leaving the source at ``time_on_source``.

    Parameters
    ----------
    pre_sample, beamstop:
        Both monitors, in both models.
    positions:
        The position of every component, in m, by name, as :meth:`.Instrument.positions`.
    time_on_source:
        When every sub-pulse leaves the source, for the flight to the pre-sample monitor alone.
    model:
        Which map of the comparisons: ``"first"``, McStas, or ``"second"``, tof.
    min_fraction:
        The weakest patch that counts as a sub-pulse, as a fraction of the strongest.

    Returns
    -------
    :
        The intensity of each sub-pulse at the pre-sample monitor, along ``subpulse``, with the
        times at both monitors; the coordinates of :func:`.incident_flight` between the
        monitors, ``Ltotal``, ``tof``, ``incident_wavelength`` and ``incident_energy``, with the
        prefix ``from_two_monitors_``, and from the source, with the prefix
        ``from_single_monitor_``; and ``reference_energy``, from the mean
        wavelength at the pre-sample monitor.

    Raises
    ------
    ValueError
        If the monitors do not have as many sub-pulses, or one arrives at the beamstop monitor
        first.
    """
    pre_sample_patches, beamstop_patches = (
        find_patches(getattr(comparison, model).values, comparison.bins)
        for comparison in (pre_sample, beamstop)
    )
    pre_sample_position = sc.vector(positions["pre_sample_monitor"], unit="m")
    try:
        # --8<-- [start:two_flights]
        flights = subpulse_flight(
            pre_sample_patches,
            pre_sample_position,
            beamstop_patches,
            sc.vector(positions["beamstop_monitor"], unit="m"),
            min_fraction=min_fraction,
        )
        from_source = incident_flight(
            sc.vector(positions["source"], unit="m"),
            time_on_source.to(unit="ms"),
            pre_sample_position,
            flights.coords["time_at_start"],
        )
        # --8<-- [end:two_flights]
    except ValueError as error:
        raise ValueError(f"{model}: {error}") from error
    flights.coords["time_at_pre_sample_monitor"] = flights.coords.pop("time_at_start")
    flights.coords["time_at_beamstop_monitor"] = flights.coords.pop("time_at_end")
    flights.coords["reference_energy"] = tof_kernels.energy_from_wavelength(
        wavelength=flights.coords.pop("wavelength_at_start")
    )
    for name in FLIGHT_COORDS:
        flights.coords[f"{TWO_MONITORS}_{name}"] = flights.coords.pop(name)
        flights.coords[f"{SINGLE_MONITOR}_{name}"] = from_source[name]
    return flights


def print_incident_energies(flights: dict[str, sc.DataArray]) -> None:
    """Print the time of flight and the incident energy of each sub-pulse, in each model.

    Parameters
    ----------
    flights:
        The sub-pulses of :func:`subpulse_flights`, by the name of the model.
    """
    for model, f in flights.items():
        print(
            f"{model}: flight path {f.coords[f'{TWO_MONITORS}_Ltotal'].value:.4f} m between the "
            f"monitors, {f.coords[f'{SINGLE_MONITOR}_Ltotal'].value:.4f} m from the source"
        )
        print(
            "  sub-pulse  t_pre [ms]  tof [ms]  λi [Å]  Ei [meV]  ΔEi/Ei"
            "  Ei 1 mon. [meV]  ΔEi/Ei  Ei from λ [meV]"
        )
        reference = f.coords["reference_energy"].to(unit="meV").values
        for i, (t, flight, wavelength, ei, one, ref) in enumerate(
            zip(
                f.coords["time_at_pre_sample_monitor"].values,
                f.coords[f"{TWO_MONITORS}_tof"].to(unit="ms").values,
                f.coords[f"{TWO_MONITORS}_incident_wavelength"].to(unit="Å").values,
                f.coords[f"{TWO_MONITORS}_incident_energy"].to(unit="meV").values,
                f.coords[f"{SINGLE_MONITOR}_incident_energy"].to(unit="meV").values,
                reference,
            )
        ):
            print(
                f"  {i:9d}  {t:10.3f}  {flight:8.4f}  {wavelength:6.4f}  {ei:8.3f}"
                f"  {(ei - ref) / ref:+7.2%}  {one:15.3f}  {(one - ref) / ref:+7.2%}"
                f"  {ref:15.3f}"
            )


def main(
    out: Path | None = None,
    ncount: float = 1e8,
    neutrons: int = 10_000_000,
    seed: int = 1,
    transmission: float | None = None,
) -> list[Path]:
    """Run both models of :class:`TRexToBeamstopMonitor`, print Ei and save the figures.

    The figures are ``trex_determine_ei.png``, the time profiles, from
    :func:`.plot_time_profiles`; ``trex_determine_ei_values.png``, the incident energy of each
    sub-pulse, from :func:`.plot_incident_energies`; and ``trex_determine_ei_graph.svg``, the graph
    of the calculation, from :func:`.save_graph`.

    Parameters
    ----------
    out:
        Folder of the figures. If ``None``, ``trex_incident_energy`` in the results folder of
        :meth:`.Workspace.find`.
    ncount:
        Number of McStas neutrons.
    neutrons:
        Number of tof neutrons.
    seed:
        Seed of the random numbers of both models.
    transmission:
        The factor on tof at both monitors, on top of the scale of the source monitor. If
        ``None``, tuned so that tof has the total of McStas at the pre-sample monitor.

    Returns
    -------
    :
        The paths of the figures: the time profiles, the incident energies, then the graph.
    """
    instrument = define_instrument()
    name = "trex_determine_ei"

    # run McStas
    sim = Simulation(instrument)
    mcstas_run = sim.run(folder_name=name, ncount=ncount, seed=seed)

    # run tof; its counts are scaled to the total of McStas at the source monitor, in n/s, and,
    # at the pre-sample and beamstop monitors, by the transmission on top of it
    tof_model = build_tof(instrument, neutrons=neutrons, seed=seed)
    tof_run = tof_model.run()
    bins = {
        name: monitor_bins(getattr(instrument, name).bin_limits)
        for name in ("source_monitor", *MONITORS)
    }
    mcstas_maps = read_mcstas_maps(mcstas_run, bins)
    tof_maps = read_tof_maps(tof_run, bins)
    scale = compare.rescale_tof_to_mcstas(tof_maps["source_monitor"], mcstas_maps["source_monitor"])
    if transmission is None:  # tuned at the pre-sample monitor, and kept at the beamstop monitor
        first = MONITORS[0]
        transmission = compare.rescale_tof_to_mcstas(tof_maps[first], mcstas_maps[first]) / scale
    comparisons = [
        compare.mcstas_vs_tof(
            name, mcstas_maps[name], tof_maps[name], ncount, neutrons, scale, transmission
        )
        for name in MONITORS
    ]
    for c in comparisons:
        print(
            f"{c.name}: McStas {c.first.values.sum():.4g} n/s, "
            f"tof {c.second.values.sum():.4g} n/s, transmission {transmission:.4g}"
        )

    # the incident energy of each sub-pulse, from its flight between the monitors and from the
    # source to the pre-sample monitor
    positions = instrument.positions()
    # every sub-pulse leaves the source when the choppers are phased to, for the pre-sample
    # monitor alone
    time_on_source = sc.scalar(instrument.params.t0_offset, unit="s")
    pre_sample, beamstop = comparisons
    flights = {
        label: subpulse_flights(pre_sample, beamstop, positions, time_on_source, model)
        for label, model in (("McStas", "first"), ("tof", "second"))
    }
    print_incident_energies(flights)

    conditions = instrument.running_conditions()
    folder = cli.output_folder(out, "trex_incident_energy")
    profiles = plot_time_profiles(comparisons, f"T-Rex direct beam\n{conditions}", MONITOR_LABELS)
    energies = plot_incident_energies(
        flights,
        f"T-Rex incident energy per sub-pulse\n{conditions}",
        methods=METHODS,
        time="time_at_pre_sample_monitor",
        time_label="Time at the pre-sample monitor [ms]",
    )
    files = [
        cli.save_figure(profiles, folder / f"{name}.png"),
        cli.save_figure(energies, folder / f"{name}_values.png"),
        save_graph(folder / f"{name}_graph"),
    ]
    print(f"Saved {files[-1]}")
    return files


if __name__ == "__main__":
    main(out=cli.parse_out(__doc__))
