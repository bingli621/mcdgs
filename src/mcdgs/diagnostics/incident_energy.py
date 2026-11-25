"""The incident energy of each sub-pulse, drawn, for any instrument.

:func:`plot_time_profiles` draws the intensity of monitors against time, in McStas and tof,
:func:`plot_incident_energies` the incident energy of each sub-pulse by each method, and
:func:`save_graph` the graph of the calculation, ``EI_GRAPH``. The incident energy diagnostic of
each instrument finds its sub-pulses and their incident energy.
"""

from collections.abc import Mapping, Sequence
from pathlib import Path

import matplotlib.pyplot as plt
import scipp as sc
from matplotlib.figure import Figure

from mcdgs.diagnostics import compare
from mcdgs.reduction.incident_energy import EI_GRAPH


def plot_time_profiles(
    comparisons: Sequence[compare.MonitorComparison],
    title: str,
    labels: Mapping[str, str] | None = None,
) -> Figure:
    """Draw the intensity of each monitor against time, summed over wavelength, in both models.

    Each monitor has its own color; McStas is a solid line, tof a dashed one.

    Parameters
    ----------
    comparisons:
        One per monitor, in beam order.
    title:
        The title of the figure.
    labels:
        The label of each monitor by its name, such as ``"Pre-sample monitor"``. A monitor
        without one is labeled with its name.

    Returns
    -------
    :
        The figure.
    """
    fig, ax = plt.subplots(figsize=(11, 6), layout="constrained")
    for color, c in zip(("tab:blue", "tab:orange"), comparisons):
        times = sc.midpoints(c.bins["toa"]).to(unit="ms").values
        label = (labels or {}).get(c.name, c.name)
        ax.plot(times, c.first.values.sum(axis=1), color=color, lw=1.2, label=f"{label}, McStas")
        ax.plot(
            times, c.second.values.sum(axis=1), color=color, lw=1.2, ls="--", label=f"{label}, tof"
        )
    # the monitors have time bins of the same width
    ylabel = compare.intensity_label(comparisons[0].first.sum("wavelength"))
    ax.set(title=title, xlabel="Time [ms]", ylabel=ylabel)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.4)
    ax.legend()
    return fig


def plot_incident_energies(
    flights: Mapping[str, sc.DataArray],
    title: str,
    *,
    methods: Mapping[str, tuple[str, str]],
    time: str,
    time_label: str,
) -> Figure:
    """Draw the incident energy of each sub-pulse against its time at a monitor.

    The top panel has the incident energy of each model by each method, and that from the mean
    wavelength of McStas; the bottom panel the relative difference of each from the mean
    wavelength of its own model. The color is the method; McStas has filled circles, tof
    larger hollow squares.

    Parameters
    ----------
    flights:
        The sub-pulses of each model, by its name, ``"McStas"`` and ``"tof"``, along
        ``subpulse``, with ``<method>_incident_energy`` of each method, ``time`` and
        ``reference_energy``, the incident energy from the mean wavelength.
    title:
        The title of the figure.
    methods:
        The color and the label of each method, by its name, such as ``"from_two_monitors"``.
    time:
        The coordinate of the time of each sub-pulse, on the horizontal axis.
    time_label:
        The label of that axis.

    Returns
    -------
    :
        The figure.
    """
    fig, (energy_ax, difference_ax) = plt.subplots(
        2, 1, figsize=(11, 8), sharex=True, height_ratios=(2, 1), layout="constrained"
    )
    fig.suptitle(title)
    for model, f in flights.items():
        marker, size, hollow = ("s", 9, True) if model == "tof" else ("o", 6, False)
        times = f.coords[time].to(unit="ms").values
        reference = f.coords["reference_energy"].to(unit="meV").values
        for method, (color, method_label) in methods.items():
            energy = f.coords[f"{method}_incident_energy"].to(unit="meV").values
            style = {
                "marker": marker,
                "markersize": size,
                "ls": "none",
                "color": color,
                "mfc": "none" if hollow else color,
                "label": f"{model}, {method_label}",
            }
            energy_ax.plot(times, energy, **style)
            difference_ax.plot(times, 100 * (energy - reference) / reference, **style)
    mcstas = flights["McStas"]
    energy_ax.plot(
        mcstas.coords[time].to(unit="ms").values,
        mcstas.coords["reference_energy"].to(unit="meV").values,
        marker="x",
        ls="none",
        color="black",
        label="McStas, from the mean wavelength",
    )
    energy_ax.set(ylabel="Incident energy [meV]")
    difference_ax.axhline(0, color="0.5", lw=0.8)
    difference_ax.set(
        xlabel=time_label,
        ylabel="Difference from λ [%]",
    )
    energy_ax.legend()
    for ax in (energy_ax, difference_ax):
        ax.grid(True, alpha=0.4)
    return fig


def save_graph(path: Path) -> Path:
    """Draw ``EI_GRAPH`` with ``sc.show_graph`` and save it as an SVG file.

    The coordinates are boxes and the kernels gray labels between them, with an arrow from each
    coordinate to the kernel that takes it. It needs graphviz, which the pixi feature ``mcstas``
    installs.

    Parameters
    ----------
    path:
        The file, without the extension, which graphviz adds.

    Returns
    -------
    :
        The path of the SVG file.
    """
    graph = sc.show_graph(EI_GRAPH, simplified=False)
    return Path(graph.render(path, format="svg", cleanup=True))
