"""The elastic resolution of T-Rex at a short wavelength, from vanadium: ``lambda0`` of 1 Å.

:class:`TRexElasticResolution` is :class:`.TRex` with the vanadium cylinder of :func:`.vanadium`
on ``sample_pos``, with the beamstop and the banana detector. Its control parameters are those
of :class:`.TRex` at ``LAMBDA0``, so the choppers turn and the beam on the sample is a train of
sub-pulses, each with its own incident energy: unlike the white beam of
:mod:`.diagnostics.trex.solid_angles`.
Vanadium scatters incoherently and elastically, in every direction, so every neutron it scatters
has an energy transfer of zero: the width of its elastic line is the energy resolution of the
instrument, and whatever lies away from it comes from the reduction or the instrument, not the
sample. The cylinder scatters only once, into the window of
:func:`.solid_angles.detector_window`, which covers the banana detector.

It runs in each mode of ``MODES``: ``"FF"``, the wide slits of the pulse-shaping and
monochromatic choppers, for high flux, ``"RR"``, the narrow ones, for high resolution, and
between them ``"RF"``, the narrow slits of the pulse-shaping choppers with the wide ones of the
monochromatic choppers, and ``"FR"``, the other way around, all with the same ``lambda0`` and
``rrm``, so the same sub-pulses. :func:`main` runs each mode in two
parts, with :class:`.Simulation`, as :mod:`.diagnostics.trex.diffraction` does: from the source
to ``sample_pos``, which dumps the beam there into a beam file, then from that beam through the
vanadium to the detector, into the folders of :func:`run_names`.

The reduction is that of :mod:`.diagnostics.trex.sample_runs`, as for diffraction:

1. :func:`.sample_runs.subpulse_timing` finds the sub-pulses at the pre-sample and beamstop
   monitors of the first part, with their incident energy and the time they reach the sample.
2. :func:`.sample_runs.read_detector` reads the neutrons of the banana detector of the second
   part, by pixel,
   and :func:`.scattered_pixels` leaves out the pixels within ``MIN_TWO_THETA`` of the beam: the
   direct beam that passes the beamstop.
3. :func:`.normalized_events` keeps the neutrons inside the arrival window of a
   sub-pulse, gives them their energy transfer, and divides each by the solid angle of its pixel
   and by the flux of its sub-pulse.
4. :func:`.energy_spectra` sums them into a histogram of the energy transfer of each sub-pulse,
   averaged over the solid angle of the pixels: per incident neutron, per steradian and per meV.

:func:`.fit_elastic_lines` fits a Gaussian to the elastic line of each sub-pulse, against the
final energy Ef, and gives its FWHM, ``FWHM_PER_SIGMA`` times its standard deviation σ.
:func:`.plot_spectra` draws the spectra with their errors and their Gaussians, against E.
:func:`.q_e_maps` histograms the neutrons of each sub-pulse in |Q| and E as well, and
:func:`.plot_q_e_maps` draws the maps of all of them in one figure. :func:`plot_resolutions` and
:func:`resolution_table` compare the FWHM of the modes.

:func:`main` is :func:`simulate`, which runs McStas, then :func:`analyze`, which reduces and
draws. Run it to do both::

    pixi run python -c "from mcdgs.diagnostics.trex.elastic_resolution_incoherent_short import main; main()"

Running the module only reduces the runs again, without McStas, with :func:`analyze`::

    pixi run python -m mcdgs.diagnostics.trex.elastic_resolution_incoherent_short
"""

from pathlib import Path
from typing import NamedTuple

import matplotlib.pyplot as plt
import numpy as np
import scipp as sc
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from mcdgs.core.components import Incoherent
from mcdgs.diagnostics import cli
from mcdgs.diagnostics.elastic_resolution import (
    FIT_FRACTION,
    FONT_SIZE,
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
from mcdgs.diagnostics.trex.solid_angles import detector_window
from mcdgs.instruments.trex.choppers import Mode
from mcdgs.instruments.trex.detectors import BANANA_RADIUS, BANANA_THETA_LIMITS
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.mcstas import Simulation, Workspace
from mcdgs.reduction.inelastic import normalized_events
from mcdgs.samples import vanadium

NAME = "trex_elastic_resolution_incoherent_short"
LAMBDA0 = 1.0  # Å, the central wavelength
MIN_TWO_THETA = 5.0  # deg: the pixels closer to the beam see the direct beam past the beamstop
MODES: tuple[Mode, ...] = ("FF", "FR", "RF", "RR")  # coarse to fine; FF, first, is the reference
# C0 to C3 from the finest resolution to the coarsest
MODE_COLORS: dict[Mode, str] = {"RR": "C0", "RF": "C1", "FR": "C2", "FF": "C3"}
MODEL_POINTS = 60  # the incident energies of the analytic model, across those of the sub-pulses
TERMS = ("pulse_shaping", "monochromatic", "sample", "total")  # of the analytic model
TERM_LABELS = {
    "pulse_shaping": "Pulse-shaping choppers",
    "monochromatic": "Monochromatic choppers",
    "sample": "Sample",
    "total": "All three",
}
TERM_COLORS = {"pulse_shaping": "C4", "monochromatic": "C5", "sample": "C7", "total": "k"}
# the calculations of the width, with their label and line style: eq. (3) of Lechner (2000) with
# the opening times of the slits, the same with those of the beam, and the model
CALCULATIONS = {
    "literature": ("Lechner eq. (3): slit openings, quadrature", "--"),
    "beam": ("Beam openings, quadrature", ":"),
    "model": ("Model: beam openings, convolution", "-"),
}
RATIO_LEGEND_ROOM = 1.4  # the y axis of the ratios of plot_resolutions, past the highest
MODE_LABELS: dict[Mode, str] = {
    "FF": "FF, high flux",
    "FR": "FR, narrow monochromatic slits",
    "RF": "RF, narrow pulse-shaping slits",
    "RR": "RR, high resolution",
}


# --8<-- [start:instrument]
def sample() -> Incoherent:
    """Return the vanadium cylinder of :func:`.vanadium`, scattering once into the detector.

    It scatters into the window of :func:`.solid_angles.detector_window` only.
    """
    return vanadium(name="sample", order=1, **detector_window())


class TRexElasticResolution(TRex):
    """T-Rex with a vanadium cylinder, from the source to the banana detector.

    The cylinder is that of :func:`sample`.

    Parameters
    ----------
    params:
        Control parameters.
    source_masked:
        As in :class:`.TRex`. True by default: the masked source of
        :func:`.sources.masked_source` spends the rays where the choppers let neutrons through.
    """

    name = "trex_elastic_resolution"  # the McStas instrument name

    def __init__(self, params: Params = DEFAULT_PARAMS, source_masked: bool = True) -> None:
        super().__init__(sample=sample(), params=params, source_masked=source_masked)


def define_instrument(mode: Mode = "FF") -> TRexElasticResolution:
    """Return the instrument of this diagnostic: :class:`TRexElasticResolution` at ``LAMBDA0``.

    Parameters
    ----------
    mode:
        The slits of the pulse-shaping and monochromatic choppers: ``"FF"``, the default, for
        high flux, ``"RR"`` for high resolution, or ``"RF"`` or ``"FR"`` between.

    Examples
    --------
    >>> define_instrument().params.lambda0
    1.0
    """
    params = Params(mode=mode, lambda0=LAMBDA0)  # the defaults otherwise, in the mode
    return TRexElasticResolution(params=params)


# --8<-- [end:instrument]


class Reduction(NamedTuple):
    """The reduction of a run of :func:`simulate`, from :func:`reduce_runs`."""

    pulses: sc.DataArray  # the sub-pulses of :func:`.sample_runs.subpulse_timing`
    spectra: sc.DataArray  # of :func:`.energy_spectra`
    maps: sc.DataArray  # of :func:`.q_e_maps`


def run_names(mode: Mode) -> tuple[str, str]:
    """Return the run folders of a mode: the part to ``sample_pos``, and the part from it.

    Examples
    --------
    >>> run_names("RR")[0]
    'trex_elastic_resolution_incoherent_short_rr_run_to_sample'
    """
    return f"{NAME}_{mode.lower()}_run_to_sample", f"{NAME}_{mode.lower()}_run_from_sample"


def simulate(
    ncount: float = 1e9, seed: int = 1, modes: tuple[Mode, ...] = MODES
) -> dict[Mode, tuple[Path, Path]]:
    """Run :class:`TRexElasticResolution` in McStas in each mode, in two parts split at ``sample_pos``.

    In each mode, the first part runs from the source to ``sample_pos``, with the monitors of the
    incident energy and the beam file; the second from that beam through the vanadium to the
    detector, both into the folders of :func:`run_names`, in the runs folder of
    :meth:`.Workspace.find`. Each replaces the folder of an earlier run.

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
        energy_spectra(events, pulses, pixels.sum()),
        q_e_maps(events, pulses, two_theta(binned), pixels),
    )


def resolution_table(lines: dict[Mode, sc.Dataset], pulses: dict[Mode, sc.DataArray]) -> str:
    """Return the FWHM and the flux of each sub-pulse in every mode, as a Markdown table.

    The n-th sub-pulse of one mode is the n-th of another: the modes share ``lambda0`` and
    ``rrm``, so the sub-pulses have the same Ei.

    Parameters
    ----------
    lines:
        The fits of :func:`.fit_elastic_lines`, by mode, with the same sub-pulses.
    pulses:
        The sub-pulses of :func:`.sample_runs.subpulse_timing`, by mode, with their flux at the
        pre-sample monitor.

    Returns
    -------
    :
        One row per sub-pulse: its Ei in the first mode, the FWHM, FWHM / Ei and flux of each
        mode, and the FWHM and the flux of each mode after the first over those of the first.
    """
    modes = list(lines)
    first = modes[0]
    header = ["Sub-pulse", f"Ei {first}"]
    for mode in modes:
        header += [f"FWHM {mode}", f"FWHM / Ei {mode}", f"Flux {mode}"]
    for mode in modes[1:]:
        header += [f"FWHM {mode} / {first}", f"Flux {mode} / {first}"]
    rows = []
    for i, ei in enumerate(lines[first].coords["incident_energy"].to(unit="meV").values):
        row = [f"{i}", f"{ei:.2f} meV"]
        for mode in modes:
            fwhm = lines[mode]["fwhm"]
            row += [
                f"{fwhm.values[i]:.3f} ± {np.sqrt(fwhm.variances[i]):.3f} meV",
                f"{100 * fwhm.values[i] / ei:.2f} %",
                f"{pulses[mode].values[i]:.2e} n/s",
            ]
        for mode in modes[1:]:
            width = lines[mode]["fwhm"].values[i] / lines[first]["fwhm"].values[i]
            flux = pulses[mode].values[i] / pulses[first].values[i]
            row += [f"{width:.2f}", f"{flux:.2f}"]
        rows.append(row)
    return cli.markdown_table(header, rows)


def model_geometry(instrument: TRexElasticResolution) -> model.ResolutionGeometry:
    """Return the geometry of the analytic model of :mod:`.elastic_resolution_model`.

    It is that of the instrument, with the banana detector at ``BANANA_RADIUS``, the cylinder of
    :func:`sample`, and the scattering angles of the banana that :func:`.scattered_pixels` keeps:
    past ``MIN_TWO_THETA`` on either side of the beam.

    Parameters
    ----------
    instrument:
        The instrument, as :func:`define_instrument` returns it in a mode.

    Returns
    -------
    :
        The geometry.
    """
    start, stop = BANANA_THETA_LIMITS
    two_theta_ranges = ((start, -MIN_TWO_THETA), (MIN_TWO_THETA, stop))
    return model.geometry_of(instrument, BANANA_RADIUS, float(sample().radius), two_theta_ranges)


def model_geometries(modes: tuple[Mode, ...] = MODES) -> dict[Mode, model.ResolutionGeometry]:
    """Return the geometry of :func:`model_geometry` of :func:`define_instrument` in each mode.

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


def model_resolutions(
    lines: dict[Mode, sc.Dataset],
    geometries: dict[Mode, model.ResolutionGeometry],
    points: int | None = MODEL_POINTS,
) -> dict[Mode, sc.Dataset]:
    """Compute the analytic FWHM of :mod:`.elastic_resolution_model` across the Ei of each mode.

    Parameters
    ----------
    lines:
        The fits of :func:`.fit_elastic_lines`, by mode: the model spans their Ei.
    geometries:
        The geometries of :func:`model_geometry`, by mode, those of ``lines`` among them.
    points:
        The incident energies of the model, evenly spaced. If ``None``, the Ei of the sub-pulses.

    Returns
    -------
    :
        The widths of :func:`.elastic_resolution_model.elastic_resolution` in ``geometries``,
        by mode.
    """
    models = {}
    for mode, fits in lines.items():
        ei = fits.coords["incident_energy"].to(unit="meV")
        if points is not None:
            ei = sc.linspace("incident_energy", ei.min(), ei.max(), points)
        models[mode] = model.elastic_resolution(geometries[mode], ei, FIT_FRACTION)
    return models


def model_table(lines: dict[Mode, sc.Dataset], models: dict[Mode, sc.Dataset]) -> str:
    """Return the FWHM of the simulated lines and of the analytic model, as a Markdown table.

    Parameters
    ----------
    lines:
        The fits of :func:`.fit_elastic_lines`, by mode, with the same sub-pulses.
    models:
        The widths of :func:`model_resolutions` at the Ei of the sub-pulses, by mode.

    Returns
    -------
    :
        One row per sub-pulse: its Ei in the first mode, then in each mode the FWHM of McStas, of
        the model, and the model over McStas.
    """
    modes = list(lines)
    first = modes[0]
    header = ["Sub-pulse", f"Ei {first}"]
    for mode in modes:
        header += [f"McStas {mode}", f"Model {mode}", f"Model / McStas {mode}"]
    rows = []
    for i, ei in enumerate(lines[first].coords["incident_energy"].to(unit="meV").values):
        row = [f"{i}", f"{ei:.2f} meV"]
        for mode in modes:
            fwhm = lines[mode]["fwhm"]
            width = models[mode]["fwhm"].values[i]
            row += [
                f"{fwhm.values[i]:.3f} ± {np.sqrt(fwhm.variances[i]):.3f} meV",
                f"{width:.3f} meV",
                f"{width / fwhm.values[i]:.3f}",
            ]
        rows.append(row)
    return cli.markdown_table(header, rows)


def model_terms(
    geometries: dict[Mode, model.ResolutionGeometry], calculation: str = "model"
) -> tuple[float, dict[Mode, sc.Dataset]]:
    """Compute the FWHM of each term at the Ei of ``lambda0``, in every mode.

    Parameters
    ----------
    geometries:
        The geometries of :func:`model_geometry`, by mode, with the same ``lambda0``.
    calculation:
        ``"model"``, :func:`.elastic_resolution_model.elastic_resolution`, or ``"literature"``,
        :func:`.elastic_resolution_model.literature_resolution` with the openings of the slits.

    Returns
    -------
    :
        The Ei of ``lambda0``, in meV, and the widths there, by mode.

    Raises
    ------
    ValueError
        If ``calculation`` is neither ``"model"`` nor ``"literature"``.
    """
    if calculation not in ("model", "literature"):
        raise ValueError(f"calculation {calculation!r} is not 'model' or 'literature'")
    ei = model.energy_of(next(iter(geometries.values())).lambda0)
    energies = sc.array(dims=["incident_energy"], values=[ei], unit="meV")
    if calculation == "literature":
        return ei, {
            mode: model.literature_resolution(geometry, energies)
            for mode, geometry in geometries.items()
        }
    return ei, {
        mode: model.elastic_resolution(geometry, energies, FIT_FRACTION)
        for mode, geometry in geometries.items()
    }


def terms_table(ei: float, terms: dict[Mode, sc.Dataset]) -> str:
    """Return the FWHM of each term of the analytic model at one Ei, as a Markdown table.

    Parameters
    ----------
    ei:
        The incident energy, in meV.
    terms:
        The widths of :func:`model_terms`, by mode.

    Returns
    -------
    :
        One row per mode: the FWHM of the pulse-shaping choppers, of the monochromatic choppers
        and of the sample, each alone, and of the three together, in meV and as a fraction of Ei.
    """
    header = ["Mode", "Pulse-shaping", "Monochromatic", "Sample", "All three", "All three / Ei"]
    rows = []
    for mode, widths in terms.items():
        # the dataset names the width of the three together fwhm
        values = [widths["fwhm" if name == "total" else name].values[0] for name in TERMS]
        rows.append(
            [mode, *(f"{value:.3f} meV" for value in values), f"{100 * values[-1] / ei:.2f} %"]
        )
    return cli.markdown_table(header, rows)


def plot_model_terms(
    ei: float, geometries: dict[Mode, model.ResolutionGeometry], title: str
) -> Figure:
    """Draw the distribution of δE of each term of the analytic model, and of the three together.

    One panel per mode, at one Ei, each distribution over its own highest value, so that their
    widths compare. The legend gives the FWHM of each, as :func:`.elastic_resolution_model.fwhm_of`
    fits it.

    Parameters
    ----------
    ei:
        The incident energy, in meV.
    geometries:
        The geometries of :func:`model_geometry`, by mode, one panel each.
    title:
        The title of the figure.

    Returns
    -------
    :
        The figure.
    """
    modes = list(geometries)
    columns = 2
    rows = -(-len(modes) // columns)
    with plt.rc_context({"font.size": FONT_SIZE}):
        fig, axes = plt.subplots(
            rows,
            columns,
            figsize=(16, 5.5 * rows),
            sharex=True,
            squeeze=False,
            layout="constrained",
        )
        fig.suptitle(title)
        for ax, mode in zip(axes.flat, modes):
            shapes = model.line_shapes(geometries[mode], ei)
            for name in TERMS:
                energies, density = shapes[name]
                fwhm = model.fwhm_of(energies, density, FIT_FRACTION)
                ax.plot(
                    energies,
                    density / density.max(),
                    color=TERM_COLORS[name],
                    lw=2.5 if name == "total" else 1.5,
                    label=f"{TERM_LABELS[name]}: FWHM {fwhm:.3f} meV",
                )
            ax.set(title=MODE_LABELS.get(mode, mode), ylabel="Over its highest value")
            ax.set_ylim(0, 1.45)  # the legend above the lines
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize="small", loc="upper left")
        for ax in axes.flat[len(modes) :]:
            ax.set_visible(False)
        for ax in axes[-1]:
            ax.set(xlabel="δE, the error of the energy transfer [meV]")
            ax.set_xlim(-2.5, 2.5)
    return fig


def literature_resolutions(
    lines: dict[Mode, sc.Dataset],
    geometries: dict[Mode, model.ResolutionGeometry],
    points: int | None = MODEL_POINTS,
    openings: str = "slit",
) -> dict[Mode, sc.Dataset]:
    """Compute the FWHM of the formula of the literature across the Ei of each mode.

    Parameters
    ----------
    lines:
        The fits of :func:`.fit_elastic_lines`, by mode: the curves span their Ei.
    geometries:
        As :func:`model_resolutions`.
    points:
        As :func:`model_resolutions`.
    openings:
        As :func:`.elastic_resolution_model.literature_resolution`.

    Returns
    -------
    :
        The widths of :func:`.elastic_resolution_model.literature_resolution` in
        ``geometries``, by mode.
    """
    curves = {}
    for mode, fits in lines.items():
        ei = fits.coords["incident_energy"].to(unit="meV")
        if points is not None:
            ei = sc.linspace("incident_energy", ei.min(), ei.max(), points)
        curves[mode] = model.literature_resolution(geometries[mode], ei, openings)
    return curves


def calculations(
    lines: dict[Mode, sc.Dataset],
    geometries: dict[Mode, model.ResolutionGeometry],
    points: int | None,
) -> dict[str, dict[Mode, sc.Dataset]]:
    """Compute the FWHM of the three calculations of ``CALCULATIONS`` across the Ei of each mode.

    Parameters
    ----------
    lines:
        The fits of :func:`.fit_elastic_lines`, by mode.
    geometries:
        As :func:`model_resolutions`.
    points:
        As :func:`model_resolutions`.

    Returns
    -------
    :
        By calculation, as ``CALCULATIONS`` names them, the widths by mode.
    """
    return {
        "literature": literature_resolutions(lines, geometries, points, "slit"),
        "beam": literature_resolutions(lines, geometries, points, "beam"),
        "model": model_resolutions(lines, geometries, points),
    }


def plot_calculations(
    lines: dict[Mode, sc.Dataset],
    curves: dict[str, dict[Mode, sc.Dataset]],
    at_subpulses: dict[str, dict[Mode, sc.Dataset]],
    title: str,
) -> Figure:
    """Draw the FWHM of McStas and of the three calculations of ``CALCULATIONS``, in every mode.

    One column per mode. The top row draws the FWHM of the simulated lines, points with the error
    of their fit, and of each calculation, a line; the bottom row each calculation over McStas at
    the Ei of the sub-pulses, with the error of the fit.

    Parameters
    ----------
    lines:
        The fits of :func:`.fit_elastic_lines`, by mode.
    curves:
        The widths of :func:`calculations` across the Ei, by calculation and mode.
    at_subpulses:
        The widths of :func:`calculations` at the Ei of the sub-pulses, by calculation and mode.
    title:
        The title of the figure.

    Returns
    -------
    :
        The figure.
    """
    modes = list(lines)
    with plt.rc_context({"font.size": FONT_SIZE}):
        fig, axes = plt.subplots(
            2,
            len(modes),
            figsize=(5 * len(modes), 10),
            sharex=True,
            sharey="row",
            squeeze=False,
            layout="constrained",
        )
        fig.suptitle(title)
        for column, mode in enumerate(modes):
            top, bottom = axes[0, column], axes[1, column]
            color = MODE_COLORS.get(mode, "k")
            ei = lines[mode].coords["incident_energy"].to(unit="meV").values
            fwhm = lines[mode]["fwhm"]
            error = np.sqrt(fwhm.variances)
            top.errorbar(
                ei, fwhm.values, yerr=error, fmt="o", ms=6, capsize=3, color=color, label="McStas"
            )
            for name, (label, style) in CALCULATIONS.items():
                curve = curves[name][mode]
                top.plot(
                    curve.coords["incident_energy"].values,
                    curve["fwhm"].values,
                    color="k",
                    linestyle=style,
                    lw=1.5,
                    label=label,
                )
                ratio = at_subpulses[name][mode]["fwhm"].values / fwhm.values
                bottom.errorbar(
                    ei,
                    ratio,
                    yerr=ratio * error / fwhm.values,
                    fmt="o",
                    ms=5,
                    capsize=2,
                    color="k",
                    markerfacecolor="none" if name != "model" else "k",
                    linestyle=style,
                    lw=1,
                    label=label,
                )
            bottom.axhline(1, color=color, lw=1.5)
            top.set(title=MODE_LABELS.get(mode, mode))
            bottom.set(xlabel="Ei [meV]")
            for ax in (top, bottom):
                ax.grid(True, alpha=0.3)
            top.set_ylim(bottom=0)
        axes[0, 0].set(ylabel="FWHM [meV]")
        axes[1, 0].set(ylabel="Calculation / McStas")
        axes[0, 0].legend(fontsize="small", loc="upper left")
    return fig


def calculations_table(
    lines: dict[Mode, sc.Dataset], at_subpulses: dict[str, dict[Mode, sc.Dataset]]
) -> str:
    """Return each calculation of ``CALCULATIONS`` over McStas at each sub-pulse, as Markdown.

    Parameters
    ----------
    lines:
        The fits of :func:`.fit_elastic_lines`, by mode, with the same sub-pulses.
    at_subpulses:
        The widths of :func:`calculations` at the Ei of the sub-pulses, by calculation and mode.

    Returns
    -------
    :
        One row per sub-pulse: its Ei in the first mode, then in each mode each calculation over
        the FWHM of McStas.
    """
    modes = list(lines)
    first = modes[0]
    short = {"literature": "Lechner", "beam": "Beam, quadrature", "model": "Model"}
    header = ["Sub-pulse", f"Ei {first}"]
    for mode in modes:
        header += [f"{short[name]} / McStas {mode}" for name in CALCULATIONS]
    rows = []
    for i, ei in enumerate(lines[first].coords["incident_energy"].to(unit="meV").values):
        row = [f"{i}", f"{ei:.2f} meV"]
        for mode in modes:
            mcstas = lines[mode]["fwhm"].values[i]
            row += [
                f"{at_subpulses[name][mode]['fwhm'].values[i] / mcstas:.3f}"
                for name in CALCULATIONS
            ]
        rows.append(row)
    return cli.markdown_table(header, rows)


def lechner_table(
    pulses: dict[Mode, sc.DataArray], geometries: dict[Mode, model.ResolutionGeometry]
) -> str:
    """Return the opening times, their ratio and the intensity rule of Lechner (2000), as Markdown.

    At the sub-pulse nearest the Ei of ``lambda0``: the FWHM ``τ1`` and ``τ2`` of the openings of
    the pulse-shaping and monochromatic pairs, :func:`.elastic_resolution_model.opening_of`, their
    ratio ``ρ`` against the optimum ``ρ_opt`` of eq. (4), and the intensity over that of the first
    mode three ways: ``τ1 τ2`` as eq. (1) has it, the product of the areas of the two openings, and
    the flux of McStas at the pre-sample monitor.

    Parameters
    ----------
    pulses:
        The sub-pulses of :func:`.sample_runs.subpulse_timing`, by mode, with their flux.
    geometries:
        The geometries of :func:`model_geometry`, by mode, those of ``pulses`` among them.

    Returns
    -------
    :
        One row per mode.
    """
    modes = list(pulses)
    openings = {}
    for mode in modes:
        geometry = geometries[mode]
        openings[mode] = [model.opening_of(getattr(geometry, name)) for name in model.PAIRS]
    ei = model.energy_of(geometries[modes[0]].lambda0)
    first = pulses[modes[0]]
    index = int(np.argmin(np.abs(first.coords["incident_energy"].to(unit="meV").values - ei)))
    (tau1_first, area1_first), (tau2_first, area2_first) = openings[modes[0]]
    header = [
        "Mode",
        "τ1",
        "τ2",
        "ρ = τ1 / τ2",
        "ρ_opt",
        f"τ1 τ2 / {modes[0]}",
        f"Areas / {modes[0]}",
        f"McStas flux / {modes[0]}",
    ]
    rows = []
    for mode in modes:
        (tau1, area1), (tau2, area2) = openings[mode]
        flux = pulses[mode].values[index] / first.values[index]
        rows.append(
            [
                mode,
                f"{tau1 * 1e6:.0f} µs",
                f"{tau2 * 1e6:.1f} µs",
                f"{tau1 / tau2:.1f}",
                f"{model.optimal_ratio(geometries[mode]):.1f}",
                f"{tau1 * tau2 / (tau1_first * tau2_first):.2f}",
                f"{area1 * area2 / (area1_first * area2_first):.2f}",
                f"{flux:.2f}",
            ]
        )
    return cli.markdown_table(header, rows)


def plot_resolutions(
    lines: dict[Mode, sc.Dataset],
    pulses: dict[Mode, sc.DataArray],
    title: str,
    models: dict[Mode, sc.Dataset] | None = None,
) -> Figure:
    """Draw the FWHM of the elastic line and the flux of each sub-pulse against Ei, in every mode.

    Four panels on one Ei axis, each with one linear y axis: the FWHM, the FWHM as a fraction of Ei,
    the flux at the pre-sample monitor, and the FWHM and the flux of each mode after the first over
    those of the first, the price of the resolution. Each sub-pulse is a point, not joined, in the
    color of its mode, ``MODE_COLORS``; in the last panel an open circle is the FWHM and an open
    square the flux. The FWHM has the error of the fit of :func:`.fit_elastic_lines`; the flux has
    none: :func:`.find_patches` sums the intensity of the monitor without its errors, which are
    small from 1e9 rays. The n-th sub-pulse of one mode is the n-th of another, as in
    :func:`resolution_table`. The first panel also draws the FWHM of ``models``, the analytic
    model of :mod:`.elastic_resolution_model`, as a line in the color of its mode.

    Parameters
    ----------
    lines:
        The fits of :func:`.fit_elastic_lines`, by mode, with the same sub-pulses.
    pulses:
        The sub-pulses of :func:`.sample_runs.subpulse_timing`, by mode, with their flux.
    title:
        The title of the figure.
    models:
        The widths of :func:`model_resolutions`, by mode. If ``None``, no lines.

    Returns
    -------
    :
        The figure.
    """
    modes = list(lines)
    first = modes[0]
    with plt.rc_context({"font.size": FONT_SIZE}):
        fig, axes = plt.subplots(2, 2, figsize=(16, 11), sharex=True, layout="constrained")
        fig.suptitle(title)
        (absolute, relative), (flux_panel, ratios) = axes
        for mode in modes:
            color = MODE_COLORS.get(mode, "k")
            ei = lines[mode].coords["incident_energy"].to(unit="meV").values
            fwhm = lines[mode]["fwhm"]
            error = np.sqrt(fwhm.variances)
            style = {"color": color, "fmt": "o", "ms": 7, "capsize": 3}
            label = MODE_LABELS.get(mode, mode)
            absolute.errorbar(ei, fwhm.values, yerr=error, label=label, **style)
            if models is not None and mode in models:
                curve = models[mode]
                absolute.plot(
                    curve.coords["incident_energy"].values,
                    curve["fwhm"].values,
                    color=color,
                    lw=1.5,
                )
            relative.errorbar(ei, 100 * fwhm.values / ei, yerr=100 * error / ei, **style)
            flux_ei = pulses[mode].coords["incident_energy"].to(unit="meV").values
            flux = pulses[mode].to(unit="1/s").values / 1e6
            flux_panel.plot(flux_ei, flux, linestyle="none", marker="o", ms=7, color=color)
        ei = lines[first].coords["incident_energy"].to(unit="meV").values
        for mode in modes[1:]:
            color = MODE_COLORS.get(mode, "k")
            width = lines[mode]["fwhm"].data / lines[first]["fwhm"].data  # with the fit errors
            flux = pulses[mode].values / pulses[first].values
            ratios.errorbar(
                ei,
                width.values,
                yerr=np.sqrt(width.variances),
                color=color,
                fmt="o",
                ms=7,
                markerfacecolor="none",
                capsize=3,
                label=f"FWHM {mode} / {first}",
            )
            ratios.plot(
                ei,
                flux,
                linestyle="none",
                marker="s",
                ms=7,
                markerfacecolor="none",
                color=color,
                label=f"Flux {mode} / {first}",
            )
        absolute.set(ylabel="FWHM [meV]", title="The width of the elastic line")
        relative.set(ylabel="FWHM / Ei [%]", title="As a fraction of Ei")
        flux_panel.set(
            xlabel="Ei [meV]",
            ylabel="Neutrons through the whole monitor [$10^6$ n/s]",
            title="The flux at the pre-sample monitor",
        )
        ratios.set(
            xlabel="Ei [meV]",
            ylabel=f"Ratio to {first}",
            title=f"The price of the resolution, over {first}",
        )
        for ax in axes.flat:
            ax.grid(True, alpha=0.3)
            ax.set_ylim(bottom=0)
        ratios.set_ylim(top=RATIO_LEGEND_ROOM * ratios.get_ylim()[1])  # the legend above the data
        handles, labels = absolute.get_legend_handles_labels()
        if models:
            handles.append(Line2D([], [], color="gray", lw=1.5))
            labels.append("Lines: the analytic model")
        absolute.legend(handles, labels, loc="upper left")
        handles, labels = ratios.get_legend_handles_labels()  # the FWHM in one column, the flux
        order = sorted(range(len(labels)), key=lambda i: not labels[i].startswith("FWHM"))
        ratios.legend(
            [handles[i] for i in order],
            [labels[i] for i in order],
            ncols=2,
            fontsize="small",
            loc="upper center",
        )
    return fig


def analyze(out: Path | None = None, modes: tuple[Mode, ...] = MODES) -> list[Path]:
    """Reduce the runs of every mode, save their spectra and fits, and compare the modes.

    It runs no McStas, so it can run again on the same runs; :func:`main` runs them first. For
    each mode, ``m`` its name in lower case, the files are
    ``NAME_m.png``, from :func:`.plot_spectra`, ``NAME_m_lines.txt``, the table of
    :func:`.lines_table`, and ``NAME_m_sqe.png``, from :func:`.plot_q_e_maps`, with ``NAME``
    ``trex_elastic_resolution_incoherent_short``; with more than one mode, ``NAME_modes.png``,
    from :func:`plot_resolutions`, and ``NAME_modes.txt``, the table of :func:`resolution_table`.
    The analytic model adds ``NAME_model.txt``, the table of :func:`model_table`,
    ``NAME_terms.txt``, that of :func:`terms_table`, and ``NAME_terms.png``, from
    :func:`plot_model_terms`; and, against the literature, ``NAME_literature_terms.txt``, the
    table of :func:`terms_table` for the formula of the literature, ``NAME_calculations.png``,
    from :func:`plot_calculations`, ``NAME_calculations.txt``, the table of
    :func:`calculations_table`, and ``NAME_lechner.txt``, that of :func:`lechner_table`.

    Parameters
    ----------
    out:
        Folder of the files. If ``None``, ``trex_elastic_resolution_incoherent_short`` in the
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
