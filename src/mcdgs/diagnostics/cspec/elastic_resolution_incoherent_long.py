"""The elastic resolution of CSPEC at a long wavelength, from vanadium: ``lambda0`` of 8 Å.

:mod:`.diagnostics.cspec.elastic_resolution_incoherent_short` at ``LAMBDA0`` instead of the
3 Å of :class:`.Params`: the same instrument, :class:`.CSPECElasticResolution` with the vanadium
cylinder and the detector panels, the same two runs split at ``sample_pos``, and the same
reduction, fits and figures, into the folders and files of ``NAME``. Only ``lambda0`` changes:
the choppers are phased to it, and the source simulates its band.

:func:`main` is :func:`simulate`, which runs McStas, then :func:`analyze`, which reduces and
draws. Run it to do both::

    pixi run python -c "from mcdgs.diagnostics.cspec.elastic_resolution_incoherent_long import main; main()"

Running the module only reduces the runs again, without McStas, with :func:`analyze`::

    pixi run python -m mcdgs.diagnostics.cspec.elastic_resolution_incoherent_long
"""

from pathlib import Path

from mcdgs.diagnostics import cli
from mcdgs.diagnostics.cspec import sample_runs as runs
from mcdgs.diagnostics.cspec.elastic_resolution_incoherent_short import (
    CSPECElasticResolution,
    Reduction,
    plot_resolution,
)
from mcdgs.diagnostics.elastic_resolution import (
    energy_spectra,
    fit_elastic_lines,
    lines_table,
    plot_q_e_maps,
    plot_spectra,
    q_e_maps,
    two_theta,
)
from mcdgs.instruments.cspec.model import Params
from mcdgs.mcstas import Simulation, Workspace
from mcdgs.reduction.inelastic import normalized_events

NAME = "cspec_elastic_resolution_incoherent_long"
LAMBDA0 = 8.0  # Å, the central wavelength


# --8<-- [start:instrument]
def define_instrument() -> CSPECElasticResolution:
    """Return the instrument of this diagnostic: :class:`.CSPECElasticResolution` at ``LAMBDA0``.

    Examples
    --------
    >>> define_instrument().params.lambda0
    8.0
    """
    params = Params(lambda0=LAMBDA0)  # the defaults otherwise
    return CSPECElasticResolution(params=params)


# --8<-- [end:instrument]


def run_names() -> tuple[str, str]:
    """Return the run folders: the part to ``sample_pos``, and the part from it.

    Examples
    --------
    >>> run_names()[0]
    'cspec_elastic_resolution_incoherent_long_run_to_sample'
    """
    return f"{NAME}_run_to_sample", f"{NAME}_run_from_sample"


def simulate(ncount: float = 1e9, seed: int = 1) -> tuple[Path, Path]:
    """Run the instrument of :func:`define_instrument` in McStas, in two parts at ``sample_pos``.

    As :func:`.elastic_resolution_incoherent_short.simulate`, into the folders of
    :func:`run_names`.

    Parameters
    ----------
    ncount:
        Number of McStas rays of the first part, from the masked source. The second part runs
        every neutron of the beam file.
    seed:
        Seed of the random numbers of the first part; the second part takes ``seed + 1``.

    Returns
    -------
    :
        The folders of the first and the second part.
    """
    to_name, from_name = run_names()
    sim = Simulation(define_instrument())
    to_sample = sim.run_to_sample(folder_name=to_name, ncount=ncount, seed=seed)
    # sample_pos -> detector, from the beam file, in the frame of sample_pos: the beam frame
    from_sample = sim.run_from_sample(folder_name=from_name, seed=seed + 1)
    return to_sample.folder, from_sample.folder


def reduce_runs(
    to_sample: str | Path | None = None, from_sample: str | Path | None = None
) -> Reduction:
    """Reduce the two parts of a McStas run of :func:`simulate`, without running McStas.

    As :func:`.elastic_resolution_incoherent_short.reduce_runs`.

    Parameters
    ----------
    to_sample, from_sample:
        The folders of the first and the second part. If ``None``, those of :func:`run_names`
        in the runs folder of :meth:`.Workspace.find`, where :func:`simulate` writes.

    Returns
    -------
    :
        The sub-pulses of :func:`.sample_runs.subpulse_timing`, the spectra of
        :func:`.energy_spectra` and the maps of :func:`.q_e_maps`.
    """
    workspace = Workspace.find()
    to_name, from_name = run_names()
    to_sample = Path(to_sample or workspace.run_folder(to_name))
    from_sample = Path(from_sample or workspace.run_folder(from_name))
    pulses = runs.subpulse_timing_of_run(to_sample, define_instrument())  # as simulate runs it
    binned = runs.read_detector(from_sample)
    pixels = runs.pixel_solid_angle(binned)
    events = normalized_events(binned, pulses, pixels)
    return Reduction(
        pulses,
        energy_spectra(events, pulses, pixels.sum()),
        q_e_maps(events, pulses, two_theta(binned), pixels),
    )


def analyze(out: Path | None = None) -> list[Path]:
    """Reduce the runs, and save their spectra, their fits and their S(Q, E) maps.

    It runs no McStas, so it can run again on the same runs; :func:`main` runs them first. The
    files are ``cspec_elastic_resolution_incoherent_long.png``, from :func:`.plot_spectra`,
    ``cspec_elastic_resolution_incoherent_long_lines.txt``, the table of :func:`.lines_table`,
    ``cspec_elastic_resolution_incoherent_long_sqe.png``, from :func:`.plot_q_e_maps`, and
    ``cspec_elastic_resolution_incoherent_long_resolution.png``, from :func:`.plot_resolution`.

    Parameters
    ----------
    out:
        Folder of the files. If ``None``, ``cspec_elastic_resolution_incoherent_long`` in the
        results folder of :meth:`.Workspace.find`.

    Returns
    -------
    :
        The paths of the files.
    """
    folder = cli.output_folder(out, NAME)
    pulses, spectra, maps = reduce_runs()
    conditions = define_instrument().running_conditions()
    lines = fit_elastic_lines(spectra)
    table = lines_table(lines)
    print(table, end="")
    files = [cli.save_text(table, folder / f"{NAME}_lines.txt")]
    fig = plot_spectra(
        spectra,
        lines,
        f"CSPEC, vanadium: the energy transfer of {pulses.sizes['subpulse']} sub-pulses\n"
        f"{conditions}",
    )
    files.append(cli.save_figure(fig, folder / f"{NAME}.png"))
    fig = plot_q_e_maps(maps, f"CSPEC, vanadium: S(Q, E) of each sub-pulse\n{conditions}")
    files.append(cli.save_figure(fig, folder / f"{NAME}_sqe.png"))
    fig = plot_resolution(
        lines,
        pulses,
        f"CSPEC, vanadium: the resolution and the flux of each sub-pulse\n{conditions}",
    )
    files.append(cli.save_figure(fig, folder / f"{NAME}_resolution.png"))
    return files


def main(out: Path | None = None, ncount: float = 1e9, seed: int = 1) -> list[Path]:
    """Run McStas with :func:`simulate`, then reduce and draw with :func:`analyze`.

    Parameters
    ----------
    out:
        Folder of the files, as :func:`analyze`.
    ncount, seed:
        As :func:`simulate`.

    Returns
    -------
    :
        The paths of the files.
    """
    simulate(ncount=ncount, seed=seed)
    return analyze(out)


if __name__ == "__main__":
    analyze(out=cli.parse_out(__doc__))  # McStas runs only with main
