"""The beam at the entrance of guide 1 of T-Rex, from McStas and from tof.

:class:`TRexSource` is the source and a time and wavelength monitor at the entrance of
guide 1, the first optics of T-Rex, 1.895927 m from the source with the opening of that guide.
Nothing stands between the two, so McStas and tof should see the same pulse there: only their
models of the ESS moderator differ.

Counts and intensity. Each simulated neutron stands for some of the real flux, its weight:
intensity = counts × weight, in n/s. The counts are how many simulated neutrons land in a bin,
and depend only on how many are simulated; the intensity is the physics, and should not.

- McStas sets the weights itself, so that the rays add up to the flux of the moderator. With
  ten times the ``ncount``, a bin has ten times the rays, each of a tenth of the weight: the
  intensity stays, only its noise goes down.
- tof's pulse is a distribution with no flux: each neutron has weight 1, so its intensity is
  its count, and grows with the number of neutrons drawn. So tof takes its weight from McStas:
  ``scale``, the total of McStas over that of tof, is the n/s each tof neutron stands for.
  Scaled, the level is that of McStas, and the shape in time and wavelength that of each
  model.

Error bars. A bin of N neutrons has the Poisson error √N, a relative error of 1/√N. In n/s it
is √N × weight. A hundred times the neutrons gives ten times the error in counts, but a tenth
of the relative error, and a tenth of the error in n/s, as the weight is a hundredth. So the
error bars shrink with more neutrons, in both models, while the intensity stays.

:func:`main` runs both models, reads the intensity of each on the bins of the monitor, scales
tof to McStas, and :func:`plot_comparison` draws them side by side. Run the module to do all
of it::

    pixi run python -m mcdgs.diagnostics.trex.source
"""

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import mcstasscript as ms
import numpy as np
import scipp as sc
from matplotlib.figure import Figure
from numpy.typing import NDArray

from mcdgs.core.base import Instrument
from mcdgs.core.components import ESS_butterfly, Monitor_nD
from mcdgs.core.params import InstrParam, param
from mcdgs.instruments.trex.model import DEFAULT_PARAMS as TREX_DEFAULTS
from mcdgs.instruments.trex.sources import source
from mcdgs.mcstas.simulation import Simulation
from mcdgs.tof.builder import build_tof
from mcdgs.tof.reader import read_tof_monitor

# --8<-- [start:trex_source]
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
# --8<-- [end:trex_source]


def plot_comparison(
    mcstas_intensity: NDArray,
    mcstas_error: NDArray,
    tof_intensity: NDArray,
    tof_error: NDArray,
    bins: Mapping[str, sc.Variable],
    ncounts: float,
    neutrons: int,
    scale: float,
) -> Figure:
    """Draw both models at the entrance of guide 1, in n/s.

    The top row has the time and wavelength map of each, on one color scale; the bottom row
    the wavelength spectrum and the time profile of both, summed over the other axis, with
    their error bars.

    Parameters
    ----------
    mcstas_intensity:
        The intensity of McStas per bin, in n/s, with the shape (time bins, wavelength bins).
    mcstas_error:
        Its error per bin, with the same shape. The bins are independent, so summed over an
        axis their errors add in quadrature.
    tof_intensity:
        The intensity of tof per bin, in n/s, with the same shape: its counts scaled to the
        total of McStas.
    tof_error:
        Its error per bin, with the same shape: the Poisson error of the counts, √N, scaled
        by the same factor. Summed over an axis, they add in quadrature too.
    bins:
        The bin edges, ``"toa"`` and ``"wavelength"``, evenly spaced, in units of time and of
        length.
    ncounts:
        Number of McStas neutrons, for the legend.
    neutrons:
        Number of tof neutrons, for the legend.
    scale:
        The factor that scales tof's counts to n/s, for the legend.

    Returns
    -------
    :
        The figure.
    """
    times = bins["toa"].to(unit="ms").values
    wavelengths = bins["wavelength"].to(unit="Å").values
    time_bin, wavelength_bin = times[1] - times[0], wavelengths[1] - wavelengths[0]
    vmax = max(mcstas_intensity.max(), tof_intensity.max())  # one color scale for both maps

    fig, ((mcstas_ax, tof_ax), (spectrum_ax, profile_ax)) = plt.subplots(
        2, 2, figsize=(12, 9), layout="constrained"
    )
    fig.suptitle("T-Rex Flux at entrance of guide 1")

    # the time and wavelength maps
    mcstas_ax.pcolormesh(times, wavelengths, mcstas_intensity.T, cmap="turbo", vmin=0, vmax=vmax)
    mesh = tof_ax.pcolormesh(times, wavelengths, tof_intensity.T, cmap="turbo", vmin=0, vmax=vmax)
    mcstas_ax.set(title="McStas", xlabel="Time [ms]", ylabel="Wavelength [Å]")
    tof_ax.set(title="tof", xlabel="Time [ms]", ylabel="Wavelength [Å]")
    fig.colorbar(
        mesh,
        ax=[mcstas_ax, tof_ax],
        label=f"Intensity [n/s per {time_bin:.2g} ms and {wavelength_bin:.2g} Å]",
    )

    time_centers = sc.midpoints(bins["toa"]).to(unit="ms").values
    wavelength_centers = sc.midpoints(bins["wavelength"]).to(unit="Å").values
    errorbar_style = {"fmt": ".", "markersize": 3, "elinewidth": 0.8}
    mcstas_label = f"McStas, ncount = {int(ncounts):,}"
    tof_label = f"tof, {neutrons:,} neutrons × {scale:.3g} n/s each"

    # the wavelength spectrum, summed over time
    spectrum_ax.errorbar(
        wavelength_centers,
        mcstas_intensity.sum(axis=0),
        yerr=np.sqrt((mcstas_error**2).sum(axis=0)),
        label=mcstas_label,
        **errorbar_style,
    )
    spectrum_ax.errorbar(
        wavelength_centers,
        tof_intensity.sum(axis=0),
        yerr=np.sqrt((tof_error**2).sum(axis=0)),
        label=tof_label,
        **errorbar_style,
    )
    spectrum_ax.set(
        xlabel="Wavelength [Å]", ylabel=f"Intensity [n/s per {wavelength_bin:.2g} Å]", ylim=0
    )
    spectrum_ax.legend()

    # the time profile, summed over wavelength
    profile_ax.errorbar(
        time_centers,
        mcstas_intensity.sum(axis=1),
        yerr=np.sqrt((mcstas_error**2).sum(axis=1)),
        label=mcstas_label,
        **errorbar_style,
    )
    profile_ax.errorbar(
        time_centers,
        tof_intensity.sum(axis=1),
        yerr=np.sqrt((tof_error**2).sum(axis=1)),
        label=tof_label,
        **errorbar_style,
    )
    profile_ax.set(xlabel="Time [ms]", ylabel=f"Intensity [n/s per {time_bin:.2g} ms]", ylim=0)
    profile_ax.legend()
    return fig


def _mean(intensity: NDArray, centers: NDArray, axis: int) -> float:
    """Return the mean of a bin center, weighted by the intensity summed over the other axis."""
    return float(np.average(centers, weights=intensity.sum(axis=1 - axis)))


def main(
    ncount: float = 1e7,
    neutrons: int = 1_000_000,
    seed: int = 1,
) -> Path:
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

    Returns
    -------
    :
        The path of the figure.
    """
    instrument = TRexSource()

    # --8<-- [start:run_mcstas]
    # run McStas and read its monitor, in n/s per bin
    folder = Simulation(instrument).run("trex_source", ncount=ncount, seed=seed).parent
    monitor = ms.name_search("source_monitor", ms.load_data(str(folder)))
    # McStasScript has the map as (wavelength, time), the y axis first
    mcstas_intensity = np.asarray(monitor.Intensity, dtype=float).T
    mcstas_error = np.asarray(monitor.Error, dtype=float).T
    # --8<-- [end:run_mcstas]

    # --8<-- [start:run_tof]
    # run tof and read its monitor, in counts per bin
    result = build_tof(instrument, neutrons=neutrons, seed=seed).run()
    bins = {  # the bins of the McStas monitor, in the order of its map: time, wavelength
        "toa": sc.array(dims=["toa"], values=TIME_EDGES, unit="s"),
        "wavelength": sc.array(dims=["wavelength"], values=WAVELENGTH_EDGES, unit="Å"),
    }
    tof_counts = read_tof_monitor(result["source_monitor"], bins)
    # --8<-- [end:run_tof]

    # print the mean time and wavelength of each model
    time_centers = sc.midpoints(bins["toa"]).to(unit="ms").values
    wavelength_centers = sc.midpoints(bins["wavelength"]).values  # Å
    for label, intensity in (("McStas", mcstas_intensity), ("tof", tof_counts)):
        print(
            f"{label:>6}: mean time {_mean(intensity, time_centers, 0):.3f} ms, "
            f"mean wavelength {_mean(intensity, wavelength_centers, 1):.3f} Å"
        )

    # scale tof to the total of McStas, in n/s, and print that total
    scale = mcstas_intensity.sum() / tof_counts.sum()
    tof_intensity = tof_counts * scale
    tof_error = np.sqrt(tof_counts) * scale  # Poisson, scaled as the counts
    print(
        f"McStas total {mcstas_intensity.sum():.4g} n/s, which tof's {tof_counts.sum():.0f} are scaled to"
    )

    figure = folder / "trex_source.png"
    fig = plot_comparison(
        mcstas_intensity, mcstas_error, tof_intensity, tof_error, bins, ncount, neutrons, scale
    )
    fig.savefig(figure, dpi=150)
    plt.close(fig)
    print(f"Saved {figure}")
    return figure


if __name__ == "__main__":
    main(ncount=1e6, neutrons=10_000)
