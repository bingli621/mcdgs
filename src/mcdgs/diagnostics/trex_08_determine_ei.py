"""The incident energy of T-Rex, from the direct beam at the pre-sample and beamstop monitors.

:class:`TRexToBeamstopMonitor` is :class:`.TRex` with the La2Ni7 crystal of :func:`.la2ni7` on
``sample_pos``, without the beamstop and the banana detector: the beamstop monitor is the last
monitor that counts. Its control parameters are those of :class:`.TRex`.

The beamstop monitor is ``BEAMSTOP_DISTANCE`` after the sample in space, but before it in beam
order: it counts the direct beam as if there were no sample. So both monitors see the same
sub-pulses, the beamstop monitor later by the flight over the distance between them.

As in :mod:`.trex_07_neutrons_on_sample`, tof's counts are scaled at the source monitor, to the
total of McStas there. A transmission on top of that scale stands for what the guides take out of
the beam, which tof does not model: by default, it is tuned so that tof has the total of McStas at
the pre-sample monitor, and the same factor holds at the beamstop monitor.

The incident energy of each sub-pulse follows from its flight between the monitors, and also
from its flight from the source to the pre-sample monitor alone, leaving the source at
``t0_offset``, the time the choppers are phased to.
:func:`subpulse_flights` finds the sub-pulses at each monitor with :func:`.find_patches` and pairs
them; ``EI_GRAPH``, built from the kernels of scippneutron, turns the start and end of a straight
flight, and the times there, into the flight path, the time of flight, the incident wavelength
and the incident energy. The mean wavelength of each sub-pulse at the pre-sample monitor, which the
monitors record, gives the incident energy it should have.

:func:`plot_time_profiles` draws the intensity of both monitors against time, and
:func:`plot_incident_energies` the incident energy of each sub-pulse. Run the module to do it::

    pixi run python -m mcdgs.diagnostics.trex_08_determine_ei
"""

from pathlib import Path
from typing import cast

import matplotlib.pyplot as plt
import mcstasscript as ms
import scipp as sc
from matplotlib.figure import Figure
from numpy.typing import NDArray
from scipp.coords.graph import GraphDict
from scippneutron.conversion import tof as tof_kernels
from scippneutron.conversion.graph import beamline

from mcdgs.core.components import Single_crystal
from mcdgs.diagnostics.cli import output_folder, parse_out
from mcdgs.diagnostics.compare import MonitorComparison, read_scaled_comparisons
from mcdgs.diagnostics.trex_07_neutrons_on_sample import find_patches
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.mcstas.simulation import Simulation
from mcdgs.samples import la2ni7
from mcdgs.tof.builder import build_tof
from mcdgs.tof.reader import monitor_bins

MONITORS = ("pre_sample_monitor", "beamstop_monitor")  # in beam order
MONITOR_LABELS = {
    "pre_sample_monitor": "Pre-sample monitor",
    "beamstop_monitor": "Beamstop monitor",
}
MIN_FRACTION = 0.02  # a patch weaker than this fraction of the strongest is not a sub-pulse


def time_of_flight(time_at_start: sc.Variable, time_at_end: sc.Variable) -> sc.Variable:
    """Return the time of flight from the start to the end of a straight flight."""
    return time_at_end - time_at_start


# a straight flight from a start to an end, as scippneutron's beam from a source to a position,
# without scattering: from the pre-sample to the beamstop monitor, or from the source to the
# pre-sample monitor
EI_GRAPH: GraphDict = {
    **beamline.beamline(scatter=False),  # Ltotal, from source_position to position
    "source_position": "start_position",
    "position": "end_position",
    "tof": time_of_flight,
    "incident_wavelength": tof_kernels.wavelength_from_tof,  # from tof and Ltotal
    "incident_energy": tof_kernels.incident_energy_from_wavelength,
}


def incident_flight(
    start_position: sc.Variable,
    time_at_start: sc.Variable,
    end_position: sc.Variable,
    time_at_end: sc.Variable,
) -> dict[str, sc.Variable]:
    """Compute the incident energy of a straight flight with ``EI_GRAPH``.

    Parameters
    ----------
    start_position, end_position:
        Where the flight starts and ends, as vectors.
    time_at_start, time_at_end:
        When it passes each, along ``subpulse`` or, for the start, one time for all.

    Returns
    -------
    :
        ``Ltotal``, ``tof``, ``incident_wavelength`` and ``incident_energy``, by name.
    """
    flight = sc.DataArray(
        sc.zeros(sizes=time_at_end.sizes),  # only the coordinates count
        coords={
            "start_position": start_position,
            "time_at_start": time_at_start,
            "end_position": end_position,
            "time_at_end": time_at_end,
        },
    )
    names = ["Ltotal", "tof", "incident_wavelength", "incident_energy"]
    flight = cast(
        sc.DataArray, flight.transform_coords(names, graph=EI_GRAPH, keep_intermediate=False)
    )
    return {name: flight.coords[name] for name in names}


class TRexToBeamstopMonitor(TRex):
    """T-Rex with a La2Ni7 crystal, from the source to the beamstop monitor.

    The beamstop and the banana detector, which only see what the sample scatters or lets
    through, are left out.

    Parameters
    ----------
    sample:
        The crystal. If ``None``, :func:`.la2ni7`.
    params:
        Control parameters.
    """

    name = "trex_to_beamstop_monitor"  # the McStas instrument name

    def __init__(
        self, sample: Single_crystal | None = None, params: Params = DEFAULT_PARAMS
    ) -> None:
        super().__init__(sample=la2ni7(name="sample") if sample is None else sample, params=params)
        # after the sample: the beamstop monitor is the last one that counts
        del self.beamstop
        del self.banana_detector


def subpulse_flights(
    pre_sample: MonitorComparison,
    beamstop: MonitorComparison,
    positions: dict[str, list[float]],
    time_on_source: sc.Variable,
    model: str,
    min_fraction: float = MIN_FRACTION,
) -> sc.DataArray:
    """Pair the sub-pulses of one model at both monitors, and compute their incident energy.

    The sub-pulses are the patches of :func:`.find_patches` with at least ``min_fraction`` of
    the intensity of the strongest. The flight between the monitors is shorter than the time
    between sub-pulses, so the n-th sub-pulse at one is the n-th at the other.

    The incident energy comes from two flights, each through :func:`incident_flight`: from the
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
        ``"mcstas"`` or ``"tof"``.
    min_fraction:
        The weakest patch that counts as a sub-pulse, as a fraction of the strongest.

    Returns
    -------
    :
        The intensity of each sub-pulse at the pre-sample monitor, along ``subpulse``, with the
        times at both monitors; the coordinates of :func:`incident_flight` between the
        monitors, ``Ltotal``, ``tof``, ``incident_wavelength`` and ``incident_energy``, and from
        the source, with the prefix ``one_monitor_``; and ``reference_energy``, from the mean
        wavelength at the pre-sample monitor.

    Raises
    ------
    ValueError
        If the monitors do not have as many sub-pulses, or one arrives at the beamstop monitor
        first.
    """
    found = []
    for comparison in (pre_sample, beamstop):
        intensity: NDArray = getattr(comparison, f"{model}_intensity")
        patches = find_patches(intensity, comparison.bins)
        keep = patches.intensity >= min_fraction * patches.intensity.max()
        found.append((patches, keep))
    (first, keep_first), (second, keep_second) = found
    if keep_first.sum() != keep_second.sum():
        raise ValueError(
            f"{model}: {keep_first.sum()} sub-pulses at the pre-sample monitor, "
            f"{keep_second.sum()} at the beamstop monitor"
        )
    time_at_pre_sample = sc.array(dims=["subpulse"], values=first.time[keep_first], unit="ms")
    time_at_beamstop = sc.array(dims=["subpulse"], values=second.time[keep_second], unit="ms")
    pre_sample_position = sc.vector(positions["pre_sample_monitor"], unit="m")
    between_monitors = incident_flight(
        pre_sample_position,
        time_at_pre_sample,
        sc.vector(positions["beamstop_monitor"], unit="m"),
        time_at_beamstop,
    )
    from_source = incident_flight(
        sc.vector(positions["source"], unit="m"),
        time_on_source.to(unit="ms"),
        pre_sample_position,
        time_at_pre_sample,
    )
    flights = sc.DataArray(
        sc.array(dims=["subpulse"], values=first.intensity[keep_first], unit="1/s"),  # n/s
        coords={
            "time_at_pre_sample_monitor": time_at_pre_sample,
            "time_at_beamstop_monitor": time_at_beamstop,
            **between_monitors,
            **{f"one_monitor_{name}": value for name, value in from_source.items()},
            "reference_energy": tof_kernels.energy_from_wavelength(
                wavelength=sc.array(
                    dims=["subpulse"], values=first.wavelength[keep_first], unit="Å"
                )
            ),
        },
    )
    if (flights.coords["tof"] <= sc.scalar(0.0, unit="ms")).any().value:
        raise ValueError(f"{model}: a sub-pulse reaches the beamstop monitor first")
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
            f"{model}: flight path {f.coords['Ltotal'].value:.4f} m between the monitors, "
            f"{f.coords['one_monitor_Ltotal'].value:.4f} m from the source"
        )
        print(
            "  sub-pulse  t_pre [ms]  tof [ms]  λi [Å]  Ei [meV]  ΔEi/Ei"
            "  Ei 1 mon. [meV]  ΔEi/Ei  Ei from λ [meV]"
        )
        reference = f.coords["reference_energy"].to(unit="meV").values
        for i, (t, flight, wavelength, ei, one, ref) in enumerate(
            zip(
                f.coords["time_at_pre_sample_monitor"].values,
                f.coords["tof"].to(unit="ms").values,
                f.coords["incident_wavelength"].to(unit="Å").values,
                f.coords["incident_energy"].to(unit="meV").values,
                f.coords["one_monitor_incident_energy"].to(unit="meV").values,
                reference,
            )
        ):
            print(
                f"  {i:9d}  {t:10.3f}  {flight:8.4f}  {wavelength:6.4f}  {ei:8.3f}"
                f"  {(ei - ref) / ref:+7.2%}  {one:15.3f}  {(one - ref) / ref:+7.2%}"
                f"  {ref:15.3f}"
            )


def plot_time_profiles(comparisons: list[MonitorComparison], title: str) -> Figure:
    """Draw the intensity of each monitor against time, summed over wavelength, in both models.

    Each monitor has its own color; McStas is a solid line, tof a dashed one.

    Parameters
    ----------
    comparisons:
        One per monitor, in beam order.
    title:
        The title of the figure.

    Returns
    -------
    :
        The figure.
    """
    fig, ax = plt.subplots(figsize=(11, 6), layout="constrained")
    for color, c in zip(("tab:blue", "tab:orange"), comparisons):
        times = sc.midpoints(c.bins["toa"]).to(unit="ms").values
        label = MONITOR_LABELS.get(c.name, c.name)
        ax.plot(
            times, c.mcstas_intensity.sum(axis=1), color=color, lw=1.2, label=f"{label}, McStas"
        )
        ax.plot(
            times, c.tof_intensity.sum(axis=1), color=color, lw=1.2, ls="--", label=f"{label}, tof"
        )
    # the monitors have time bins of the same width
    edges = comparisons[0].bins["toa"].to(unit="ms").values
    time_bin = edges[1] - edges[0]
    ax.set(title=title, xlabel="Time [ms]", ylabel=f"Intensity [n/s per {time_bin:.2g} ms]")
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.4)
    ax.legend()
    return fig


def plot_incident_energies(flights: dict[str, sc.DataArray], title: str) -> Figure:
    """Draw the incident energy of each sub-pulse against its time at the pre-sample monitor.

    The top panel has the incident energy of each model from the flight between the monitors
    and from the pre-sample monitor alone, and that from the mean wavelength of McStas; the
    bottom panel the relative difference of each from the mean wavelength of its own model.
    The color is the method; McStas has filled circles, tof larger hollow squares.

    Parameters
    ----------
    flights:
        The sub-pulses of :func:`subpulse_flights`, by the name of the model: ``"McStas"`` and
        ``"tof"``.
    title:
        The title of the figure.

    Returns
    -------
    :
        The figure.
    """
    fig, (energy_ax, difference_ax) = plt.subplots(
        2, 1, figsize=(11, 8), sharex=True, height_ratios=(2, 1), layout="constrained"
    )
    fig.suptitle(title)
    methods = {  # the coordinate of each method, its color and its label
        "incident_energy": ("tab:blue", "between the monitors"),
        "one_monitor_incident_energy": ("tab:orange", "pre-sample monitor alone"),
    }
    for model, f in flights.items():
        marker, size, hollow = ("s", 9, True) if model == "tof" else ("o", 6, False)
        times = f.coords["time_at_pre_sample_monitor"].to(unit="ms").values
        reference = f.coords["reference_energy"].to(unit="meV").values
        for coord, (color, method) in methods.items():
            energy = f.coords[coord].to(unit="meV").values
            style = {
                "marker": marker,
                "markersize": size,
                "ls": "none",
                "color": color,
                "mfc": "none" if hollow else color,
                "label": f"{model}, {method}",
            }
            energy_ax.plot(times, energy, **style)
            difference_ax.plot(times, 100 * (energy - reference) / reference, **style)
    mcstas = flights["McStas"]
    energy_ax.plot(
        mcstas.coords["time_at_pre_sample_monitor"].to(unit="ms").values,
        mcstas.coords["reference_energy"].to(unit="meV").values,
        marker="x",
        ls="none",
        color="black",
        label="McStas, from the mean wavelength",
    )
    energy_ax.set(ylabel="Incident energy [meV]")
    difference_ax.axhline(0, color="0.5", lw=0.8)
    difference_ax.set(
        xlabel="Time at the pre-sample monitor [ms]",
        ylabel="Difference from λ [%]",
    )
    energy_ax.legend()
    for ax in (energy_ax, difference_ax):
        ax.grid(True, alpha=0.4)
    return fig


def main(
    out: Path | None = None,
    ncount: float = 1e8,
    neutrons: int = 10_000_000,
    seed: int = 1,
    transmission: float | None = None,
) -> list[Path]:
    """Run both models of :class:`TRexToBeamstopMonitor`, print Ei and save the figures.

    The figures are ``trex_determine_ei.png``, the time profiles, from
    :func:`plot_time_profiles`, and ``trex_determine_ei_values.png``, the incident energy of each
    sub-pulse, from :func:`plot_incident_energies`.

    Parameters
    ----------
    out:
        Folder of the figures. If ``None``, the folder of the McStas run, next to its
        ``mccode.h5``.
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
        The paths of the figures: the time profiles, then the incident energies.
    """
    instrument = TRexToBeamstopMonitor()  # the default control parameters
    name = "trex_determine_ei"

    # run McStas and load its monitors, in n/s per bin
    folder = Simulation(instrument).run(name, ncount=ncount, seed=seed).parent
    mcstas_data = ms.load_data(str(folder))

    # run tof; its counts are scaled to the total of McStas at the source monitor, in n/s, and,
    # at the pre-sample and beamstop monitors, by the transmission on top of it
    tof_result = build_tof(instrument, neutrons=neutrons, seed=seed).run()
    read = read_scaled_comparisons(
        mcstas_data,
        tof_result,
        {name: monitor_bins(getattr(instrument, name)) for name in ("source_monitor", *MONITORS)},
        ncount,
        neutrons,
        transmission=transmission,
    )
    comparisons = [read[name] for name in MONITORS]
    for c in comparisons:
        print(
            f"{c.name}: McStas {c.mcstas_intensity.sum():.4g} n/s, "
            f"tof {c.tof_intensity.sum():.4g} n/s, transmission {c.transmission:.4g}"
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
        for label, model in (("McStas", "mcstas"), ("tof", "tof"))
    }
    print_incident_energies(flights)

    conditions = instrument.running_conditions()
    figures = []
    for suffix, fig in (
        ("", plot_time_profiles(comparisons, f"T-Rex direct beam\n{conditions}")),
        (
            "_values",
            plot_incident_energies(flights, f"T-Rex incident energy per sub-pulse\n{conditions}"),
        ),
    ):
        figure = output_folder(out, folder) / f"{name}{suffix}.png"
        fig.savefig(figure, dpi=150)
        plt.close(fig)
        print(f"Saved {figure}")
        figures.append(figure)
    return figures


if __name__ == "__main__":
    main(out=parse_out(__doc__))
