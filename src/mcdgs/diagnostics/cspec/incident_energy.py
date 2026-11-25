"""The incident energy of CSPEC, from the direct beam at the monitor after the monochromatic choppers.

:class:`CSPECToMcMonitor` is :class:`.CSPEC` up to ``mc_monitor``, after ``mc2``, with time bins
of ``MC_MONITOR_TIME_RESOLUTION``, finer than those of CSPEC, to time each sub-pulse. Its control
parameters are those of :class:`.CSPEC`. CSPEC has no pre-sample or beamstop monitor yet, so the
incident energy of each sub-pulse follows from its flight from the source to ``mc_monitor``
alone, leaving the source at ``t0_offset``, the time the choppers are phased to, as the
pre-sample monitor alone gives it for T-Rex.

As in :mod:`.diagnostics.cspec.beam_monitor`, tof's counts are scaled at the source monitor, to
the total of McStas there. A transmission on top of that scale stands for what the guides take
out of the beam, which tof does not model: by default, it is tuned so that tof has the total of
McStas at ``mc_monitor``.

:func:`subpulse_flights` finds the sub-pulses at ``mc_monitor`` with :func:`.find_patches`, and
:func:`.incident_flight` turns the start and end of the flight of each, and the times there, into
the flight path, the time of flight, the incident wavelength and the incident energy, through a
graph of scippneutron's kernels. The mean wavelength of each sub-pulse at ``mc_monitor``, which
the monitor records, gives the incident energy it should have.

:func:`.plot_time_profiles` draws the intensity of the monitor against time,
:func:`.plot_incident_energies` the incident energy of each sub-pulse, and :func:`.save_graph` the
graph, all of :mod:`mcdgs.diagnostics.incident_energy`. Run the module to do it::

    pixi run python -m mcdgs.diagnostics.cspec.incident_energy
"""

from pathlib import Path

import scipp as sc
from scippneutron.conversion import tof as tof_kernels

from mcdgs.diagnostics import cli, compare
from mcdgs.diagnostics.incident_energy import plot_incident_energies, plot_time_profiles, save_graph
from mcdgs.instruments.cspec.detectors import detector_panels
from mcdgs.instruments.cspec.model import CSPEC, DEFAULT_PARAMS, Params
from mcdgs.instruments.cspec.monitors import mc_monitor
from mcdgs.instruments.cspec.sources import wavelength_band
from mcdgs.mcstas import Simulation, read_mcstas_maps
from mcdgs.reduction.binning import monitor_bins
from mcdgs.reduction.incident_energy import incident_flight
from mcdgs.reduction.patches import MIN_FRACTION, find_patches
from mcdgs.tof import build_tof, read_tof_maps

MONITOR = "mc_monitor"  # the monitor the sub-pulses are timed at
MC_MONITOR_TIME_RESOLUTION = 5e-6  # s, as the pre-sample monitor of T-Rex
MONITOR_LABELS = {MONITOR: "Monitor after the monochromatic choppers"}
# the method of the incident energy, the prefix of the coordinates of its flight
SINGLE_MONITOR = "from_single_monitor"
METHODS = {SINGLE_MONITOR: ("tab:orange", "mc_monitor alone")}  # its color and its label


# --8<-- [start:instrument]
class CSPECToMcMonitor(CSPEC):
    """CSPEC from the source to the monitor after the monochromatic choppers.

    W03-09-02, ``sample_pos``, the sample monitor and the detector panels, after ``mc_monitor``,
    are left out.
    ``mc_monitor`` has time bins of ``MC_MONITOR_TIME_RESOLUTION``, finer than those of CSPEC.

    Parameters
    ----------
    params:
        Control parameters.
    source_masked:
        As in :class:`.CSPEC`. True by default: the masked source of
        :func:`.sources.masked_source` spends the rays where the choppers let neutrons through.
    """

    name = "cspec_to_mc_monitor"  # the McStas instrument name

    def __init__(self, params: Params = DEFAULT_PARAMS, source_masked: bool = True) -> None:
        super().__init__(params=params, source_masked=source_masked)
        self.mc_monitor = mc_monitor(  # in its place, with finer time bins
            lambda_limits=wavelength_band(params.lambda0),
            time_resolution=MC_MONITOR_TIME_RESOLUTION,
        )
        # after mc_monitor: it is the last monitor that counts
        del self.guide_09_02
        del self.sample_pos
        del self.sample_monitor
        for name in detector_panels():
            delattr(self, name)


def define_instrument() -> CSPECToMcMonitor:
    """Return the instrument of this diagnostic: :class:`CSPECToMcMonitor` at its defaults."""
    params = Params()  # the default control parameters
    return CSPECToMcMonitor(params=params)


# --8<-- [end:instrument]


def subpulse_flights(
    comparison: compare.MonitorComparison,
    positions: dict[str, list[float]],
    time_on_source: sc.Variable,
    model: str,
    min_fraction: float = MIN_FRACTION,
) -> sc.DataArray:
    """Find the sub-pulses of one model at ``mc_monitor``, and compute their incident energy.

    The sub-pulses are the patches of :func:`.find_patches` with at least ``min_fraction`` of
    the intensity of the strongest. The incident energy comes from the flight from the source to
    the monitor, through :func:`.incident_flight`, every sub-pulse leaving the source at
    ``time_on_source``.

    Parameters
    ----------
    comparison:
        ``mc_monitor``, in both models.
    positions:
        The position of every component, in m, by name, as :meth:`.Instrument.positions`.
    time_on_source:
        When every sub-pulse leaves the source.
    model:
        Which map of the comparison: ``"first"``, McStas, or ``"second"``, tof.
    min_fraction:
        The weakest patch that counts as a sub-pulse, as a fraction of the strongest.

    Returns
    -------
    :
        The intensity of each sub-pulse at the monitor, along ``subpulse``, with
        ``time_at_mc_monitor``; the coordinates of :func:`.incident_flight`, ``Ltotal``,
        ``tof``, ``incident_wavelength`` and ``incident_energy``, with the prefix
        ``from_single_monitor_``; and ``reference_energy``,
        from the mean wavelength at the monitor.
    """
    patches = find_patches(getattr(comparison, model).values, comparison.bins)
    subpulses = patches.strongest(min_fraction)
    time_at_monitor = sc.array(dims=["subpulse"], values=subpulses.time, unit="ms")
    # --8<-- [start:one_flight]
    flight = incident_flight(
        sc.vector(positions["source"], unit="m"),
        time_on_source.to(unit="ms"),
        sc.vector(positions[MONITOR], unit="m"),
        time_at_monitor,
    )
    # --8<-- [end:one_flight]
    wavelength = sc.array(dims=["subpulse"], values=subpulses.wavelength, unit="Å")
    return sc.DataArray(
        sc.array(dims=["subpulse"], values=subpulses.intensity, unit="1/s"),  # n/s
        coords={
            "time_at_mc_monitor": time_at_monitor,
            "reference_energy": tof_kernels.energy_from_wavelength(wavelength=wavelength),
            **{f"{SINGLE_MONITOR}_{name}": value for name, value in flight.items()},
        },
    )


def print_incident_energies(flights: dict[str, sc.DataArray]) -> None:
    """Print the time of flight and the incident energy of each sub-pulse, in each model.

    Parameters
    ----------
    flights:
        The sub-pulses of :func:`subpulse_flights`, by the name of the model.
    """
    for model, f in flights.items():
        length = f.coords[f"{SINGLE_MONITOR}_Ltotal"].value
        print(f"{model}: flight path {length:.4f} m from the source")
        print("  sub-pulse  t_mc [ms]  tof [ms]  λi [Å]  Ei [meV]  ΔEi/Ei  Ei from λ [meV]")
        reference = f.coords["reference_energy"].to(unit="meV").values
        for i, (t, flight, wavelength, ei, ref) in enumerate(
            zip(
                f.coords["time_at_mc_monitor"].values,
                f.coords[f"{SINGLE_MONITOR}_tof"].to(unit="ms").values,
                f.coords[f"{SINGLE_MONITOR}_incident_wavelength"].to(unit="Å").values,
                f.coords[f"{SINGLE_MONITOR}_incident_energy"].to(unit="meV").values,
                reference,
            )
        ):
            print(
                f"  {i:9d}  {t:9.3f}  {flight:8.4f}  {wavelength:6.4f}  {ei:8.3f}"
                f"  {(ei - ref) / ref:+7.2%}  {ref:15.3f}"
            )


def main(
    out: Path | None = None,
    ncount: float = 1e8,
    neutrons: int = 10_000_000,
    seed: int = 1,
    transmission: float | None = None,
) -> list[Path]:
    """Run both models of :class:`CSPECToMcMonitor`, print Ei and save the figures.

    The figures are ``cspec_determine_ei.png``, the time profile, from
    :func:`.plot_time_profiles`; ``cspec_determine_ei_values.png``, the incident energy of each
    sub-pulse, from :func:`.plot_incident_energies`; and ``cspec_determine_ei_graph.svg``, the
    graph of the calculation, from :func:`.save_graph`.

    Parameters
    ----------
    out:
        Folder of the figures. If ``None``, ``cspec_incident_energy`` in the results folder of
        :meth:`.Workspace.find`.
    ncount:
        Number of McStas neutrons.
    neutrons:
        Number of tof neutrons.
    seed:
        Seed of the random numbers of both models.
    transmission:
        The factor on tof at ``mc_monitor``, on top of the scale of the source monitor. If
        ``None``, tuned so that tof has the total of McStas there.

    Returns
    -------
    :
        The paths of the figures: the time profile, the incident energies, then the graph.
    """
    instrument = define_instrument()
    name = "cspec_determine_ei"

    # run McStas
    sim = Simulation(instrument)
    mcstas_run = sim.run(folder_name=name, ncount=ncount, seed=seed)

    # run tof; its counts are scaled to the total of McStas at the source monitor, in n/s, and,
    # at mc_monitor, by the transmission on top of it
    tof_model = build_tof(instrument, neutrons=neutrons, seed=seed)
    tof_run = tof_model.run()
    bins = {
        name: monitor_bins(getattr(instrument, name).bin_limits)
        for name in ("source_monitor", MONITOR)
    }
    mcstas_maps = read_mcstas_maps(mcstas_run, bins)
    tof_maps = read_tof_maps(tof_run, bins)
    scale = compare.rescale_tof_to_mcstas(tof_maps["source_monitor"], mcstas_maps["source_monitor"])
    if transmission is None:  # tuned at mc_monitor
        transmission = (
            compare.rescale_tof_to_mcstas(tof_maps[MONITOR], mcstas_maps[MONITOR]) / scale
        )
    comparison = compare.mcstas_vs_tof(
        MONITOR, mcstas_maps[MONITOR], tof_maps[MONITOR], ncount, neutrons, scale, transmission
    )
    print(
        f"{MONITOR}: McStas {comparison.first.values.sum():.4g} n/s, "
        f"tof {comparison.second.values.sum():.4g} n/s, transmission {transmission:.4g}"
    )

    # the incident energy of each sub-pulse, from its flight from the source to mc_monitor,
    # every sub-pulse leaving the source when the choppers are phased to
    positions = instrument.positions()
    time_on_source = sc.scalar(instrument.params.t0_offset, unit="s")
    flights = {
        label: subpulse_flights(comparison, positions, time_on_source, model)
        for label, model in (("McStas", "first"), ("tof", "second"))
    }
    print_incident_energies(flights)

    conditions = instrument.running_conditions()
    folder = cli.output_folder(out, "cspec_incident_energy")
    profiles = plot_time_profiles([comparison], f"CSPEC direct beam\n{conditions}", MONITOR_LABELS)
    energies = plot_incident_energies(
        flights,
        f"CSPEC incident energy per sub-pulse\n{conditions}",
        methods=METHODS,
        time="time_at_mc_monitor",
        time_label="Time at mc_monitor [ms]",
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
