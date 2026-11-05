"""The beam at the entrance of guide 1 of T-Rex, from McStas and from tof.

:class:`TRexSource` is the source and a time and wavelength monitor at the entrance of
guide 1, the first optics of T-Rex, 1.895927 m from the source with the opening of that guide.
Nothing stands between the two, so McStas and tof should see the same pulse there: only their
models of the ESS moderator differ.

:func:`read_mcstas` and :func:`read_tof` give each model's intensity on the bins of the monitor,
:func:`normalized` puts both on a total of 1, and :func:`plot_comparison` draws them side by
side. Run the module to do all of it::

    pixi run python -m mcdgs.diagnostics.trex.source
"""

from dataclasses import dataclass
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import tof
from matplotlib.figure import Figure
from numpy.typing import NDArray

from mcdgs.core.base import Instrument
from mcdgs.core.components import ESS_butterfly, Monitor_nD
from mcdgs.core.params import InstrParam, param
from mcdgs.instruments.trex.model import DEFAULT_PARAMS as TREX_DEFAULTS
from mcdgs.instruments.trex.sources import source
from mcdgs.mcstas.simulation import Simulation
from mcdgs.tof.builder import build_tof

TIME_EDGES = np.linspace(0.0, 0.015, 301)  # s: the whole pulse and the flight up to 10 Å
WAVELENGTH_EDGES = np.linspace(0.0, 10.0, 201)  # Å: the band of the source, 0.1 to 10 Å


@dataclass(frozen=True)
class TRexSourceParams:
    """Control parameters of :class:`TRexSource`, with the defaults of :class:`.TRex`.

    Frozen: change them with ``TRexSource.update_params``. Each field's unit and comment are in
    its metadata. Nothing uses them yet: they set no McStas parameter and no component.

    Parameters
    ----------
    lambda0:
        Central wavelength, in Å.
    rrm:
        Repetition rate multiplication factor.
    """

    lambda0: float = param(TREX_DEFAULTS.lambda0, "Å", "central wavelength")
    rrm: int = param(TREX_DEFAULTS.rrm, "", "repetition rate multiplication factor")


class TRexSource(Instrument):
    """The source of T-Rex and a time and wavelength monitor at the entrance of guide 1.

    The monitor has the opening of guide 1. It is where guide 1 begins, on the axis of the source,
    unrotated.

    Parameters
    ----------
    params:
        Control parameters. If ``None``, the defaults of :class:`TRexSourceParams`.
    """

    params_class = TRexSourceParams
    name = "trex_source"  # the McStas instrument name

    source: ESS_butterfly
    source_monitor: Monitor_nD

    def __init__(self, params: TRexSourceParams | None = None) -> None:
        super().__init__(params)
        self.source = source()
        self.source_monitor = Monitor_nD(
            xwidth=0.089955,  # m, cell E13, the entrance of guide 1
            yheight=0.029992,  # m, cell F13, the entrance of guide 1
            restore_neutron=1,
            options=(  # time in s, lambda in Å: one 2D map
                f"time limits=[{TIME_EDGES[0]} {TIME_EDGES[-1]}] bins={len(TIME_EDGES) - 1} "
                f"lambda limits=[{WAVELENGTH_EDGES[0]} {WAVELENGTH_EDGES[-1]}] "
                f"bins={len(WAVELENGTH_EDGES) - 1}"
            ),
            filename="source_monitor",  # the name of the monitor, in both models
            at=[0, 0, 1.895927],  # m, cell B13, the entrance of guide 1
            relative="source",
        )

    @property
    def instr_params(self) -> list[InstrParam]:
        """McStas instrument parameters: there are none."""
        return []


def read_mcstas(h5_file: Path) -> NDArray:
    """Read the intensity of the monitor of a McStas run of :class:`TRexSource`.

    Parameters
    ----------
    h5_file:
        The file of the run.

    Returns
    -------
    :
        The intensity per bin, in n/s, with the shape (time bins, wavelength bins) of
        :data:`TIME_EDGES` and :data:`WAVELENGTH_EDGES`.
    """
    with h5py.File(h5_file) as file:
        components = file["entry1/instrument/components"]
        assert isinstance(components, h5py.Group)
        monitor = next(c for name, c in components.items() if name.endswith("_source_monitor"))
        return np.asarray(monitor["output/source_monitor_t_L/data"][()], dtype=float)


def read_tof(result: tof.Result) -> NDArray:
    """Count the neutrons of the first pulse of a tof run that reach the monitor.

    Parameters
    ----------
    result:
        The result of a tof run of :class:`TRexSource`.

    Returns
    -------
    :
        The neutrons per bin, with the shape (time bins, wavelength bins) of :data:`TIME_EDGES`
        and :data:`WAVELENGTH_EDGES`.
    """
    neutrons = result["source_monitor"].data["pulse", 0]
    reached = neutrons[~neutrons.masks["blocked_by_others"]]
    times = reached.coords["toa"].to(unit="s").values
    wavelengths = reached.coords["wavelength"].to(unit="angstrom").values
    counts, _, _ = np.histogram2d(times, wavelengths, bins=[TIME_EDGES, WAVELENGTH_EDGES])
    return counts


def normalized(intensity: NDArray) -> NDArray:
    """Scale an intensity to a total of 1, so that the two models compare in shape only.

    McStas gives n/s and tof a count of neutrons, so only their shapes compare.

    Parameters
    ----------
    intensity:
        The intensity per bin.

    Returns
    -------
    :
        The intensity divided by its total.

    Examples
    --------
    >>> normalized(np.array([1.0, 3.0])).tolist()
    [0.25, 0.75]
    """
    total = intensity.sum()
    if total <= 0:
        raise ValueError("no intensity reaches the monitor")
    return intensity / total


def plot_comparison(mcstas: NDArray, tof_counts: NDArray) -> Figure:
    """Draw both models at the entrance of guide 1, each scaled to a total of 1.

    The top row has the time and wavelength map of each, on one color scale; the bottom row
    the wavelength spectrum and the time profile of both, summed over the other axis.

    Parameters
    ----------
    mcstas:
        The intensity of McStas, from :func:`read_mcstas`.
    tof_counts:
        The neutrons of tof, from :func:`read_tof`.

    Returns
    -------
    :
        The figure.
    """
    maps = {"McStas": normalized(mcstas), "tof": normalized(tof_counts)}
    times_ms = TIME_EDGES * 1e3
    vmax = max(m.max() for m in maps.values())

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout="constrained")
    for ax, (label, intensity) in zip(axes[0], maps.items()):
        mesh = ax.pcolormesh(
            times_ms, WAVELENGTH_EDGES, intensity.T, cmap="turbo", vmin=0, vmax=vmax
        )
        ax.set_title(f"{label} at the entrance of guide 1")
        ax.set_xlabel("Time [ms]")
        ax.set_ylabel("Wavelength [Å]")
    fig.colorbar(mesh, ax=axes[0], label="Fraction of the intensity per bin")

    spectrum_ax, profile_ax = axes[1]
    for label, intensity in maps.items():
        spectrum_ax.stairs(intensity.sum(axis=0), WAVELENGTH_EDGES, label=label)
        profile_ax.stairs(intensity.sum(axis=1), times_ms, label=label)
    spectrum_ax.set_xlabel("Wavelength [Å]")
    spectrum_ax.set_ylabel("Fraction of the intensity per bin")
    profile_ax.set_xlabel("Time [ms]")
    profile_ax.set_ylabel("Fraction of the intensity per bin")
    for ax in (spectrum_ax, profile_ax):
        ax.set_ylim(bottom=0)
        ax.legend()
    return fig


def _mean(intensity: NDArray, centers: NDArray, axis: int) -> float:
    """Return the mean of a bin center, weighted by the intensity summed over the other axis."""
    return float(np.average(centers, weights=intensity.sum(axis=1 - axis)))


def main(
    ncount: float = 1e7,
    neutrons: int = 1_000_000,
    seed: int = 1,
) -> None:
    """Run both models of :class:`TRexSource`, print their means and save the figure.

    The figure, ``trex_source.png``, goes to the folder of the McStas run, next to its
    ``mccode.h5``.

    Parameters
    ----------
    ncount:
        Number of McStas neutrons.
    neutrons:
        Number of tof neutrons.
    seed:
        Seed of the random numbers of both models.
    """
    instrument = TRexSource()
    h5_file = Simulation(instrument).run("trex_source", ncount=ncount, seed=seed)
    mcstas = read_mcstas(h5_file)
    tof_counts = read_tof(build_tof(instrument, neutrons=neutrons, seed=seed).run())

    time_centers = (TIME_EDGES[:-1] + TIME_EDGES[1:]) / 2 * 1e3  # ms
    wavelength_centers = (WAVELENGTH_EDGES[:-1] + WAVELENGTH_EDGES[1:]) / 2
    for label, intensity in (("McStas", mcstas), ("tof", tof_counts)):
        print(
            f"{label:>6}: mean time {_mean(intensity, time_centers, 0):.3f} ms, "
            f"mean wavelength {_mean(intensity, wavelength_centers, 1):.3f} Å"
        )

    figure = h5_file.parent / "trex_source.png"
    fig = plot_comparison(mcstas, tof_counts)
    fig.savefig(figure, dpi=150)
    plt.close(fig)
    print(f"Saved {figure}")


if __name__ == "__main__":
    main(ncount=1e7, neutrons=1_000_000)
