"""The neutrons that reach the sample of CSPEC, from ESS_butterfly and from Polygon_ESS_butterfly.

Both runs are :class:`CSPECToSampleMonitor`, CSPEC up to its sample monitor, with the same number
of rays and the same seed. The first has the source of :func:`.sources.source`,
``ESS_butterfly``. The second has that of :func:`.sources.masked_source`,
``Polygon_ESS_butterfly``, masked by the eight choppers. The sample monitor of both has
wavelength bins of ``LAMBDA_RESOLUTION`` and time bins of ``TIME_RESOLUTION``, finer than those
of CSPEC, to show the shape of each sub-pulse. :mod:`mcdgs.diagnostics.masked_source` runs and
compares them.

:func:`main` runs both and saves the figures and the table. Run the module to do it::

    pixi run python -m mcdgs.diagnostics.cspec.masked_source
"""

from pathlib import Path

from mcdgs.diagnostics import cli
from mcdgs.diagnostics.masked_source import (
    LAMBDA_RESOLUTION,
    TIME_RESOLUTION,
    find_peaks,
    plot_peaks,
    plot_sources,
    run_source,
    summary_table,
)
from mcdgs.instruments.cspec.model import CSPEC, DEFAULT_PARAMS, Params
from mcdgs.instruments.cspec.monitors import sample_monitor
from mcdgs.instruments.cspec.sources import wavelength_band


# --8<-- [start:instrument]
class CSPECToSampleMonitor(CSPEC):
    """CSPEC from the source to its sample monitor, without a sample, with finer bins there.

    The sample monitor spans the band of the source, with time and wavelength bins finer than
    those of CSPEC, to show the shape of each sub-pulse.

    Parameters
    ----------
    params:
        Control parameters.
    source_masked:
        As in :class:`.CSPEC`: False for the plain source, True for the masked one.
    lambda_resolution:
        The wavelength bin of the sample monitor, in Å.
    time_resolution:
        Its time bin, in s.
    """

    name = "cspec_to_sample_monitor"  # the McStas instrument name

    def __init__(
        self,
        params: Params = DEFAULT_PARAMS,
        source_masked: bool = False,
        lambda_resolution: float = LAMBDA_RESOLUTION,
        time_resolution: float = TIME_RESOLUTION,
    ) -> None:
        super().__init__(params=params, source_masked=source_masked)
        self.sample_monitor = sample_monitor(  # in its place, with finer bins
            lambda_limits=wavelength_band(params.lambda0),
            lambda_resolution=lambda_resolution,
            time_resolution=time_resolution,
        )


def define_instrument(
    source_masked: bool = False,
    lambda_resolution: float = LAMBDA_RESOLUTION,
    time_resolution: float = TIME_RESOLUTION,
) -> CSPECToSampleMonitor:
    """Return an instrument of this diagnostic: :class:`CSPECToSampleMonitor` with one source.

    The masked one has a name of its own, so that each is compiled in a build folder of its own.

    Parameters
    ----------
    source_masked:
        True for the masked source, False for the plain one.
    lambda_resolution, time_resolution:
        As in :class:`CSPECToSampleMonitor`.

    Returns
    -------
    :
        The instrument.
    """
    params = Params()  # the default control parameters
    instrument = CSPECToSampleMonitor(
        params=params,
        source_masked=source_masked,
        lambda_resolution=lambda_resolution,
        time_resolution=time_resolution,
    )
    if source_masked:
        instrument.name = "cspec_to_sample_monitor_masked"
    return instrument


# --8<-- [end:instrument]


def main(
    out: Path | None = None,
    ncount: float = 1e8,
    seed: int = 1,
    lambda_resolution: float = LAMBDA_RESOLUTION,
    time_resolution: float = TIME_RESOLUTION,
) -> list[Path]:
    """Run both sources of :class:`CSPECToSampleMonitor`, print the table and save the files.

    The figures are ``cspec_masked_source.png`` and ``cspec_masked_source_log.png``, on a linear
    and on a log intensity axis, from :func:`.plot_sources`, and
    ``cspec_masked_source_peaks.png``, from :func:`.plot_peaks`; ``cspec_masked_source.txt`` is
    the table of :func:`.summary_table`.

    Parameters
    ----------
    out:
        Folder of the files. If ``None``, ``cspec_masked_source`` in the results folder of
        :meth:`.Workspace.find`.
    ncount:
        Number of McStas rays of each run.
    seed:
        Seed of the random numbers of both runs.
    lambda_resolution:
        The wavelength bin of the sample monitor, in Å.
    time_resolution:
        Its time bin, in s.

    Returns
    -------
    :
        The paths of the figures, linear, log then the peaks, and of the table.
    """
    name = "cspec_masked_source"
    plain_instrument, masked_instrument = (
        define_instrument(masked, lambda_resolution, time_resolution) for masked in (False, True)
    )
    plain = run_source(plain_instrument, f"{name}_plain", ncount, seed)
    masked = run_source(masked_instrument, f"{name}_masked", ncount, seed)

    table = summary_table(plain, masked, ncount)
    print(table, end="")
    folder = cli.output_folder(out, "cspec_masked_source")

    peaks = find_peaks(masked)
    conditions = f"ncount = {int(ncount):,} each\n{plain_instrument.running_conditions()}"
    files = [
        cli.save_figure(
            plot_sources(
                (plain, masked),
                f"CSPEC flux at the sample monitor, {conditions}",
                log_y=log_y,
            ),
            folder / f"{name}{suffix}.png",
        )
        for suffix, log_y in (("", False), ("_log", True))
    ]
    fig = plot_peaks(
        (plain, masked),
        f"The {len(peaks)} sub-pulses at the sample monitor, {conditions}",
        peaks,
    )
    files.append(cli.save_figure(fig, folder / f"{name}_peaks.png"))
    files.append(cli.save_text(table, folder / f"{name}.txt"))
    return files


if __name__ == "__main__":
    main(out=cli.parse_out(__doc__))
