"""The elastic resolution of T-Rex at a middle wavelength, from vanadium: ``lambda0`` of 2 Å.

:mod:`.diagnostics.trex.elastic_resolution_incoherent_short` at ``LAMBDA0`` instead of its 1 Å:
the same instrument, :class:`.TRexElasticResolution` with the vanadium cylinder, the beamstop
and the banana detector, the same modes of ``MODES``, the same two runs split at ``sample_pos``
in each mode, and the same reduction, fits, analytic model and figures, into the folders and
files of ``NAME``. Only ``lambda0`` changes, to which the choppers are phased and whose band the
source simulates: the bender stays out, so the beam goes through the thermal guide 4a, as at
1 Å. The sub-pulses have about a quarter of the Ei of those at 1 Å, and their elastic lines are
0.25 to 3.4 meV wide, so the spectra and the S(Q, E) maps are binned in E by ``ENERGY_STEP``, a
fifth of the 0.25 meV at 1 Å.

:func:`main` is :func:`simulate`, which runs McStas, then :func:`analyze`, which reduces and
draws. Run it to do both::

    pixi run python -c "from mcdgs.diagnostics.trex.elastic_resolution_incoherent_middle import main; main()"

Running the module only reduces the runs again, without McStas, with :func:`analyze`::

    pixi run python -m mcdgs.diagnostics.trex.elastic_resolution_incoherent_middle
"""

from pathlib import Path

import scipp as sc

from mcdgs.diagnostics import cli
from mcdgs.diagnostics.elastic_resolution import (
    energy_spectra,
    fit_elastic_lines,
    lines_table,
    plot_q_e_maps,
    plot_spectra,
    q_e_maps,
    scattered_pixels,
    two_theta,
)
from mcdgs.diagnostics.trex import elastic_resolution_model as model
from mcdgs.diagnostics.trex import sample_runs as runs
from mcdgs.diagnostics.trex.elastic_resolution_incoherent_short import (
    MIN_TWO_THETA,
    MODEL_POINTS,
    MODES,
    Reduction,
    TRexElasticResolution,
    calculations,
    calculations_table,
    lechner_table,
    model_geometry,
    model_resolutions,
    model_table,
    model_terms,
    plot_calculations,
    plot_model_terms,
    plot_resolutions,
    resolution_table,
    terms_table,
)
from mcdgs.instruments.trex.choppers import Mode
from mcdgs.instruments.trex.model import Params
from mcdgs.mcstas import Simulation, Workspace
from mcdgs.reduction.inelastic import normalized_events

NAME = "trex_elastic_resolution_incoherent_middle"
LAMBDA0 = 2.0  # Å, the central wavelength
ENERGY_STEP = 0.05  # meV, the bin width in E of the spectra and maps: the lines are 0.25 to 3.4 meV


# --8<-- [start:instrument]
def define_instrument(mode: Mode = "FF") -> TRexElasticResolution:
    """Return the instrument: :class:`.TRexElasticResolution` at ``LAMBDA0``, with the bender out.

    Parameters
    ----------
    mode:
        The slits of the pulse-shaping and monochromatic choppers: ``"FF"``, the default, for
        high flux, ``"RR"`` for high resolution, or ``"RF"`` or ``"FR"`` between.

    Examples
    --------
    >>> define_instrument().params.lambda0
    2.0
    >>> define_instrument().params.bender_in
    False
    """
    # the defaults otherwise, in the mode: the bender out, as at 1 Å
    params = Params(mode=mode, lambda0=LAMBDA0)
    return TRexElasticResolution(params=params)


# --8<-- [end:instrument]


def run_names(mode: Mode) -> tuple[str, str]:
    """Return the run folders of a mode: the part to ``sample_pos``, and the part from it.

    Examples
    --------
    >>> run_names("RR")[0]
    'trex_elastic_resolution_incoherent_middle_rr_run_to_sample'
    """
    return f"{NAME}_{mode.lower()}_run_to_sample", f"{NAME}_{mode.lower()}_run_from_sample"


def model_geometries(modes: tuple[Mode, ...] = MODES) -> dict[Mode, model.ResolutionGeometry]:
    """Return the geometry of :func:`.model_geometry` of :func:`define_instrument` in each mode.

    Parameters
    ----------
    modes:
        The modes.

    Returns
    -------
    :
        The geometries, by mode.
    """
    return {mode: model_geometry(define_instrument(mode)) for mode in modes}


def simulate(
    ncount: float = 1e9, seed: int = 1, modes: tuple[Mode, ...] = MODES
) -> dict[Mode, tuple[Path, Path]]:
    """Run the instrument of :func:`define_instrument` in McStas in each mode, in two parts.

    As :func:`.elastic_resolution_incoherent_short.simulate`, into the folders of
    :func:`run_names`.

    Parameters
    ----------
    ncount:
        Number of McStas rays of the first part, from the masked source. The second part runs
        every neutron of the beam file.
    seed:
        Seed of the random numbers of the first part, in every mode; the second part takes
        ``seed + 1``.
    modes:
        The modes to run.

    Returns
    -------
    :
        The folders of the first and the second part, by mode.
    """
    folders = {}
    for mode in modes:
        to_name, from_name = run_names(mode)
        sim = Simulation(define_instrument(mode))
        to_sample = sim.run_to_sample(folder_name=to_name, ncount=ncount, seed=seed)
        # sample_pos -> detector, from the beam file, in the frame of sample_pos: the beam frame
        from_sample = sim.run_from_sample(folder_name=from_name, seed=seed + 1)
        folders[mode] = (to_sample.folder, from_sample.folder)
    return folders


def reduce_runs(
    mode: Mode = "FF", to_sample: str | Path | None = None, from_sample: str | Path | None = None
) -> Reduction:
    """Reduce the two parts of a McStas run of :func:`simulate`, without running McStas.

    As :func:`.elastic_resolution_incoherent_short.reduce_runs`, binned in E by ``ENERGY_STEP``.

    Parameters
    ----------
    mode:
        The mode of the run.
    to_sample, from_sample:
        The folders of the first and the second part. If ``None``, those of :func:`run_names`
        in the runs folder of :meth:`.Workspace.find`, where :func:`simulate` writes.

    Returns
    -------
    :
        The sub-pulses of :func:`.sample_runs.subpulse_timing`, the spectra of
        :func:`.energy_spectra` and the maps of :func:`.q_e_maps`, from the pixels of
        :func:`.scattered_pixels`.
    """
    workspace = Workspace.find()
    to_name, from_name = run_names(mode)
    to_sample = Path(to_sample or workspace.run_folder(to_name))
    from_sample = Path(from_sample or workspace.run_folder(from_name))
    # as simulate runs it
    pulses = runs.subpulse_timing_of_run(to_sample, define_instrument(mode))
    binned = scattered_pixels(runs.read_detector(from_sample), MIN_TWO_THETA)
    pixels = runs.pixel_solid_angle(binned)
    events = normalized_events(binned, pulses, pixels)
    return Reduction(
        pulses,
        energy_spectra(events, pulses, pixels.sum(), step=ENERGY_STEP),
        q_e_maps(events, pulses, two_theta(binned), pixels, energy_step=ENERGY_STEP),
    )


def analyze(out: Path | None = None, modes: tuple[Mode, ...] = MODES) -> list[Path]:
    """Reduce the runs of every mode, save their spectra and fits, and compare the modes.

    As :func:`.elastic_resolution_incoherent_short.analyze`, with the same files, ``NAME``
    ``trex_elastic_resolution_incoherent_middle``. It runs no McStas, so it can run again on the
    same runs; :func:`main` runs them first.

    Parameters
    ----------
    out:
        Folder of the files. If ``None``, ``trex_elastic_resolution_incoherent_middle`` in the
        results folder of :meth:`.Workspace.find`.
    modes:
        The modes to reduce, whose runs :func:`simulate` made.

    Returns
    -------
    :
        The paths of the files.
    """
    folder = cli.output_folder(out, NAME)
    geometries = model_geometries(modes)
    files = []
    lines: dict[Mode, sc.Dataset] = {}
    fluxes: dict[Mode, sc.DataArray] = {}
    for mode in modes:
        pulses, spectra, maps = reduce_runs(mode)
        fluxes[mode] = pulses
        conditions = define_instrument(mode).running_conditions()
        stem = f"{NAME}_{mode.lower()}"
        lines[mode] = fit_elastic_lines(spectra)
        table = lines_table(lines[mode])
        print(f"{mode}:\n{table}", end="")
        files.append(cli.save_text(table, folder / f"{stem}_lines.txt"))
        fig = plot_spectra(
            spectra,
            lines[mode],
            f"T-Rex, vanadium: the energy transfer of {pulses.sizes['subpulse']} sub-pulses\n"
            f"{conditions}",
        )
        files.append(cli.save_figure(fig, folder / f"{stem}.png"))
        fig = plot_q_e_maps(maps, f"T-Rex, vanadium: S(Q, E) of each sub-pulse\n{conditions}")
        files.append(cli.save_figure(fig, folder / f"{stem}_sqe.png"))
    if len(lines) > 1:
        table = resolution_table(lines, fluxes)
        print(table, end="")
        files.append(cli.save_text(table, folder / f"{NAME}_modes.txt"))
        conditions = define_instrument().running_conditions().replace(f" · {modes[0]}", "")
        fig = plot_resolutions(
            lines,
            fluxes,
            f"T-Rex, vanadium: the elastic resolution in {', '.join(lines)}\n{conditions}",
            model_resolutions(lines, geometries),
        )
        files.append(cli.save_figure(fig, folder / f"{NAME}_modes.png"))
        table = model_table(lines, model_resolutions(lines, geometries, points=None))
        print(table, end="")
        files.append(cli.save_text(table, folder / f"{NAME}_model.txt"))
        at_subpulses = calculations(lines, geometries, points=None)
        fig = plot_calculations(
            lines,
            calculations(lines, geometries, MODEL_POINTS),
            at_subpulses,
            f"T-Rex, vanadium: the elastic resolution of Lechner (2000), of the model and of "
            f"McStas\n{conditions}",
        )
        files.append(cli.save_figure(fig, folder / f"{NAME}_calculations.png"))
        table = calculations_table(lines, at_subpulses)
        print(table, end="")
        files.append(cli.save_text(table, folder / f"{NAME}_calculations.txt"))
        table = lechner_table(fluxes, geometries)
        print(table, end="")
        files.append(cli.save_text(table, folder / f"{NAME}_lechner.txt"))
    ei, terms = model_terms(geometries)
    table = terms_table(ei, terms)
    print(table, end="")
    files.append(cli.save_text(table, folder / f"{NAME}_terms.txt"))
    _, terms = model_terms(geometries, "literature")
    table = terms_table(ei, terms)
    print(table, end="")
    files.append(cli.save_text(table, folder / f"{NAME}_literature_terms.txt"))
    conditions = define_instrument().running_conditions().replace(f" · {modes[0]}", "")
    fig = plot_model_terms(
        ei,
        geometries,
        f"T-Rex, the analytic elastic resolution at Ei = {ei:.2f} meV, that of λ0\n{conditions}",
    )
    files.append(cli.save_figure(fig, folder / f"{NAME}_terms.png"))
    return files


def main(out: Path | None = None, ncount: float = 1e9, seed: int = 1) -> list[Path]:
    """Run McStas in every mode with :func:`simulate`, then reduce and draw with :func:`analyze`.

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
    analyze(out=cli.parse_out(__doc__))  # McStas runs only with main, a few minutes per mode
