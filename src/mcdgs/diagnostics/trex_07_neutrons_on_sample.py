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
``LOG_DECADES`` deep, so that the weak bins show too.

:func:`main` runs both models, reads the intensity of each on the bins of the monitors, scales
tof at the source monitor, tunes the transmission at the sample monitor and saves both figures.
Run the module to do it::

    pixi run python -m mcdgs.diagnostics.trex_07_neutrons_on_sample
"""

from pathlib import Path

import matplotlib.pyplot as plt
import mcstasscript as ms

from mcdgs.diagnostics.cli import output_folder, parse_out
from mcdgs.diagnostics.compare import plot_comparison, read_comparison
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.mcstas.simulation import Simulation
from mcdgs.tof.builder import build_tof
from mcdgs.tof.reader import monitor_bins

LOG_DECADES = 3  # the log color scale of the maps spans this many decades below the highest bin


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


def main(
    out: Path | None = None,
    ncount: float = 1e8,
    neutrons: int = 1_000_000,
    seed: int = 1,
    transmission: float | None = None,
) -> list[Path]:
    """Run both models of :class:`TRexToSample`, print the transmission and save the figures.

    The figures, from :func:`.plot_comparison`, are ``trex_neutrons_on_sample.png``, on linear
    scales, and ``trex_neutrons_on_sample_log.png``, on log scales.

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
        The paths of the figures, linear then log.
    """
    instrument = TRexToSample()  # the default control parameters
    name = "trex_neutrons_on_sample"

    # run McStas and load its monitors, in n/s per bin
    folder = Simulation(instrument).run(name, ncount=ncount, seed=seed).parent
    mcstas_data = ms.load_data(str(folder))

    # run tof, then read the monitors of both on the bins of the McStas maps: time, wavelength
    tof_result = build_tof(instrument, neutrons=neutrons, seed=seed).run()
    source_bins = monitor_bins(instrument.source_monitor)
    bins = monitor_bins(instrument.sample_monitor)
    # tof's counts are scaled to the total of McStas at the source monitor, in n/s
    at_source = read_comparison(
        "source_monitor", mcstas_data, tof_result, source_bins, ncount, neutrons
    )

    # at the sample monitor, the same scale, and the transmission on top of it
    comparison = read_comparison(
        "sample_monitor", mcstas_data, tof_result, bins, ncount, neutrons, scale=at_source.scale
    )
    if transmission is None:  # tuned: the transmission that gives tof the total of McStas
        transmission = float(comparison.mcstas_intensity.sum() / comparison.tof_intensity.sum())
    comparison = read_comparison(
        "sample_monitor",
        mcstas_data,
        tof_result,
        bins,
        ncount,
        neutrons,
        scale=at_source.scale,
        transmission=transmission,
    )

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
    return figures


if __name__ == "__main__":
    main(out=parse_out(__doc__))
