"""The elastic resolution of T-Rex, from vanadium: the intensity of each sub-pulse against E.

:class:`TRexElasticResolution` is :class:`.TRex` with the vanadium cylinder of :func:`.vanadium`
on ``sample_pos``, with the beamstop and the banana detector. Its control parameters are those
of :class:`.TRex`, so the choppers turn and the beam on the sample is a train of sub-pulses, each
with its own incident energy: unlike the white beam of :mod:`.diagnostics.trex.solid_angles`.
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
   and :func:`scattered_pixels` leaves out the pixels within ``MIN_TWO_THETA`` of the beam: the
   direct beam that passes the beamstop.
3. :func:`.normalized_events` keeps the neutrons inside the arrival window of a
   sub-pulse, gives them their energy transfer, and divides each by the solid angle of its pixel
   and by the flux of its sub-pulse.
4. :func:`energy_spectra` sums them into a histogram of the energy transfer of each sub-pulse,
   averaged over the solid angle of the pixels: per incident neutron, per steradian and per meV.

:func:`fit_elastic_lines` fits a Gaussian to the elastic line of each sub-pulse, against the
final energy Ef, and gives its FWHM, ``FWHM_PER_SIGMA`` times its standard deviation σ.
:func:`plot_spectra` draws the spectra with their errors and their Gaussians, against E.
:func:`q_e_maps` histograms the neutrons of each sub-pulse in |Q| and E as well, and
:func:`plot_q_e_maps` draws the maps of all of them in one figure. :func:`plot_resolutions` and
:func:`resolution_table` compare the FWHM of the modes.

:func:`main` is :func:`simulate`, which runs McStas, then :func:`analyze`, which reduces and
draws. Run it to do both::

    pixi run python -c "from mcdgs.diagnostics.trex.elastic_resolution_incoherent import main; main()"

Running the module only reduces the runs again, without McStas, with :func:`analyze`::

    pixi run python -m mcdgs.diagnostics.trex.elastic_resolution_incoherent
"""

from pathlib import Path
from typing import NamedTuple

import matplotlib.pyplot as plt
import numpy as np
import scipp as sc
from matplotlib.figure import Figure
from scipy.optimize import curve_fit

from mcdgs.core.components import Incoherent
from mcdgs.diagnostics import cli
from mcdgs.diagnostics.trex import sample_runs as runs
from mcdgs.diagnostics.trex.solid_angles import detector_window
from mcdgs.instruments.trex.choppers import Mode
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.mcstas import Simulation, Workspace
from mcdgs.reduction.binning import bin_edges, bins_title
from mcdgs.reduction.inelastic import DEFAULT_ENERGY_TRANSFER_RATIO, normalized_events
from mcdgs.reduction.mdnorm import momentum_transfer, q_e_coverage
from mcdgs.samples import vanadium

NAME = "trex_elastic_resolution_incoherent"
MODES: tuple[Mode, ...] = ("FF", "FR", "RF", "RR")  # coarse to fine; FF, first, is the reference
# C0 to C3 from the finest resolution to the coarsest
MODE_COLORS: dict[Mode, str] = {"RR": "C0", "RF": "C1", "FR": "C2", "FF": "C3"}
RATIO_LEGEND_ROOM = 1.4  # the y axis of the ratios of plot_resolutions, past the highest
MODE_LABELS: dict[Mode, str] = {
    "FF": "FF, high flux",
    "FR": "FR, narrow monochromatic slits",
    "RF": "RF, narrow pulse-shaping slits",
    "RR": "RR, high resolution",
}
MIN_TWO_THETA = 5.0  # deg: the pixels closer to the beam see the direct beam past the beamstop
ENERGY_STEP = 0.05  # meV, the bin width of the spectra
FWHM_PER_SIGMA = 2 * np.sqrt(2 * np.log(2))  # the full width at half maximum of a Gaussian
# the Gaussian is fitted to the bins of at least this fraction of the highest one: further out,
# a bin holds few neutrons, and the variance of their weights is no error for it
FIT_FRACTION = 0.2
SHIFT = 0.4  # each sub-pulse is lifted above the one below by this much of the highest bin
FONT_SIZE = 14  # pt, of every text of the figure
MAP_Q_STEP = 0.05  # 1/Å, the bin width in |Q| of the S(Q,E) maps; in E, ENERGY_STEP
MAP_COLUMNS = 3  # the maps of the sub-pulses, in rows of this many
# the top of the color scale of the maps, as a percentile of the bins with neutrons: a bin on the
# edge of the coverage, with little of it, can stand far above the rest
MAP_PERCENTILE = 99.9


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
    """Return the instrument of this diagnostic: :class:`TRexElasticResolution` at its defaults.

    Parameters
    ----------
    mode:
        The slits of the pulse-shaping and monochromatic choppers: ``"FF"``, the default, for
        high flux, ``"RR"`` for high resolution, or ``"RF"`` or ``"FR"`` between.
    """
    params = Params(mode=mode)  # the default control parameters, in the mode
    return TRexElasticResolution(params=params)


# --8<-- [end:instrument]


def two_theta(binned: sc.DataArray) -> sc.Variable:
    """Return the scattering angle of each pixel, from the beam along ``runs.BEAM_DIRECTION``.

    Parameters
    ----------
    binned:
        The neutrons by pixel, as :func:`.sample_runs.read_detector`, with the pixel ``position`` and the
        ``sample_position``, in the beam frame of ``sample_pos``.

    Returns
    -------
    :
        The angle between the beam and the direction from the sample to each pixel, in deg,
        along the pixel dimension.
    """
    scattered = binned.coords["position"] - binned.coords["sample_position"]
    beam = sc.vector(runs.BEAM_DIRECTION)
    cosine = sc.dot(scattered, beam) / sc.norm(scattered)
    return sc.acos(sc.values(cosine)).to(unit="deg")


def scattered_pixels(binned: sc.DataArray, min_two_theta: float = MIN_TWO_THETA) -> sc.DataArray:
    """Leave out the pixels close to the beam, which see the direct beam past the beamstop.

    Parameters
    ----------
    binned:
        The neutrons by pixel, as :func:`.sample_runs.read_detector`.
    min_two_theta:
        The smallest scattering angle kept, in deg.

    Returns
    -------
    :
        The pixels at ``two_theta`` of ``min_two_theta`` or more, from :func:`two_theta`.
    """
    (dim,) = binned.dims
    keep = two_theta(binned).values >= min_two_theta
    return binned[dim, np.flatnonzero(keep)]


def energy_spectra(
    events: sc.DataArray,
    pulses: sc.DataArray,
    solid_angle: sc.Variable,
    step: float = ENERGY_STEP,
) -> sc.DataArray:
    """Sum the neutrons of each sub-pulse into a histogram of their energy transfer.

    Each neutron of :func:`.normalized_events` times the solid angle of its pixel
    is the fraction of the incident neutrons of its sub-pulse it stands for. Summed over the
    pixels and divided by their solid angle, it is the intensity averaged over the detector.
    Each neutron of weight w counts with the variance w², so a bin has the error √Σw².

    Parameters
    ----------
    events:
        The neutrons of :func:`.normalized_events`, in 1/sr, with their
        ``subpulse``, ``energy_transfer`` and the ``solid_angle`` of their pixel.
    pulses:
        The sub-pulses of :func:`.sample_runs.subpulse_timing`.
    solid_angle:
        The solid angle of all the pixels the neutrons come from, in sr.
    step:
        The bin width of the energy transfer, in meV.

    Returns
    -------
    :
        The intensity per incident neutron, per steradian and per meV, with its variances,
        along ``subpulse`` and
        ``energy_transfer``, with the bin edges ``energy_transfer`` in meV, symmetric about 0 and
        wide enough for the arrival window of the highest Ei, and the ``incident_energy`` and
        ``incident_wavelength`` of each sub-pulse.
    """
    highest = float(pulses.coords["incident_energy"].to(unit="meV").max().value)
    limit = np.ceil(max(map(abs, DEFAULT_ENERGY_TRANSFER_RATIO)) * highest)
    edges = bin_edges("energy_transfer", -limit, limit, step, unit="meV")
    fractions = events.copy(deep=False)
    fractions.data = events.data * events.coords["solid_angle"]  # of the incident neutrons
    fractions.variances = fractions.values**2  # a neutron of weight w counts with variance w²
    fractions.coords["energy_transfer"] = events.coords["energy_transfer"].to(unit="meV")
    n = pulses.sizes["subpulse"]
    grouped = fractions.group(sc.arange("subpulse", n))
    spectra = grouped.hist(energy_transfer=edges) / solid_angle.to(unit="sr")
    spectra = spectra / sc.scalar(step, unit="meV")
    for name in ("incident_energy", "incident_wavelength"):
        spectra.coords[name] = pulses.coords[name].copy()
    return spectra


def gaussian(x: np.ndarray, height: float, center: float, sigma: float) -> np.ndarray:
    """Return a Gaussian of ``height`` at ``center``, with standard deviation ``sigma``.

    Its full width at half maximum is ``FWHM_PER_SIGMA`` times ``sigma``, about 2.355 times.

    Examples
    --------
    At half its FWHM from the center, it is at half its height:

    >>> round(float(gaussian(np.array([1.0 + FWHM_PER_SIGMA / 2]), 2.0, 1.0, 1.0)[0]), 12)
    1.0
    """
    return height * np.exp(-((x - center) ** 2) / (2 * sigma**2))


def fit_elastic_lines(spectra: sc.DataArray, fit_fraction: float = FIT_FRACTION) -> sc.Dataset:
    """Fit a Gaussian to the elastic line of each sub-pulse, against the final energy.

    The final energy of a neutron is ``Ef = Ei - E``, with the Ei of its sub-pulse. The fit
    weighs each bin by its error, and keeps only the core of the line: the bins of at least
    ``fit_fraction`` of the highest. The tails of the lines are not Gaussian, and a bin there
    holds few neutrons, some with tiny weights, whose variance is no error for the bin.

    Parameters
    ----------
    spectra:
        The spectra of :func:`energy_spectra`, with their variances.
    fit_fraction:
        The lowest bin fitted, as a fraction of the highest.

    Returns
    -------
    :
        Along ``subpulse``, each with its variance from the fit: ``center``, the Ef of the peak,
        ``sigma`` and ``fwhm``, ``FWHM_PER_SIGMA`` times ``sigma``, in meV, and ``height``; and
        ``intensity``, the spectrum summed over the energy transfer, per incident neutron and per
        steradian, and ``reduced_chi2``, the χ² of the fit per degree of freedom. The coordinates
        are those of ``spectra`` along ``subpulse``.
    """
    spectra = spectra.transpose(["subpulse", "energy_transfer"])
    edges = spectra.coords["energy_transfer"].to(unit="meV")
    transfer = sc.midpoints(edges).values
    energies = spectra.coords["incident_energy"].to(unit="meV").values
    names = ("height", "center", "sigma")
    found: dict[str, list[float]] = {f"{name}{part}": [] for name in names for part in ("", "_var")}
    chi2 = []
    for values, variances, ei in zip(spectra.values, spectra.variances, energies):
        core = (values >= fit_fraction * values.max()) & (variances > 0)
        x, y, error = ei - transfer[core], values[core], np.sqrt(variances[core])
        mean = np.average(x, weights=y)
        start = (y.max(), mean, np.sqrt(np.average((x - mean) ** 2, weights=y)))
        best, covariance = curve_fit(gaussian, x, y, p0=start, sigma=error, absolute_sigma=True)
        best[2] = abs(best[2])  # the sign of sigma does not matter
        for name, value, variance in zip(names, best, np.diag(covariance)):
            found[name].append(float(value))
            found[f"{name}_var"].append(float(variance))
        residuals = (y - gaussian(x, *best)) / error
        chi2.append(float(residuals @ residuals) / (x.size - len(names)))
    lines = sc.Dataset(
        {
            name: sc.array(
                dims=["subpulse"],
                values=found[name],
                variances=found[f"{name}_var"],
                unit=spectra.unit if name == "height" else "meV",
            )
            for name in names
        }
    )
    lines["fwhm"] = lines["sigma"] * FWHM_PER_SIGMA
    lines["intensity"] = (spectra * (edges[1:] - edges[:-1])).sum("energy_transfer").data
    lines["reduced_chi2"] = sc.array(dims=["subpulse"], values=chi2)
    for name in ("incident_energy", "incident_wavelength"):
        lines.coords[name] = spectra.coords[name].copy()
    return lines


def lines_table(lines: sc.Dataset) -> str:
    """Return the Gaussian of the elastic line of each sub-pulse as a Markdown table.

    Parameters
    ----------
    lines:
        The fits of :func:`fit_elastic_lines`.

    Returns
    -------
    :
        One row per sub-pulse: its λi and Ei, the center Ef, σ and FWHM of its Gaussian with
        their errors, the FWHM as a fraction of Ei, the χ² of the fit per degree of freedom, and
        the intensity of its spectrum.
    """

    def with_error(line: sc.DataArray, i: int, digits: int) -> str:
        return f"{line.values[i]:.{digits}f} ± {np.sqrt(line.variances[i]):.{digits}f} meV"

    rows = []
    energies = lines.coords["incident_energy"].to(unit="meV").values
    wavelengths = lines.coords["incident_wavelength"].to(unit="Å").values
    for i, (wavelength, ei) in enumerate(zip(wavelengths, energies)):
        fwhm = lines["fwhm"].values[i]
        rows.append(
            [
                f"{i}",
                f"{wavelength:.3f} Å",
                f"{ei:.2f} meV",
                with_error(lines["center"], i, 3),
                with_error(lines["sigma"], i, 3),
                with_error(lines["fwhm"], i, 3),
                f"{100 * fwhm / ei:.2f} %",
                f"{lines['reduced_chi2'].values[i]:.1f}",
                f"{lines['intensity'].values[i]:.3e}",
            ]
        )
    header = [
        "Sub-pulse",
        "λi",
        "Ei",
        "Center Ef",
        "σ",
        "FWHM",
        "FWHM / Ei",
        "χ² / dof",
        "Intensity [1/sr]",
    ]
    return cli.markdown_table(header, rows)


def _subpulse_colors(n: int) -> list[tuple[float, float, float, float]]:
    """Return a color of turbo for each of ``n`` sub-pulses, from the shortest wavelength."""
    cmap = plt.get_cmap("turbo")
    return [cmap(x) for x in np.linspace(0.05, 0.95, n)]


def plot_spectra(
    spectra: sc.DataArray,
    lines: sc.Dataset,
    title: str,
    window: tuple[float, float] = DEFAULT_ENERGY_TRANSFER_RATIO,
    shift: float = SHIFT,
) -> Figure:
    """Draw the spectrum of every sub-pulse against E, with its Gaussian.

    The scale is linear, and spans the whole energy coverage: every sub-pulse is drawn over its
    arrival window, ``window`` times its Ei, each bin as a point with its error, the empty ones
    too, and the Gaussian of :func:`fit_elastic_lines` as a line, with its FWHM. The sub-pulses
    are stacked from the highest Ei at the top to the lowest at the bottom, each lifted above
    the one below by ``shift`` of the highest bin, so that each stands on its own baseline.

    Parameters
    ----------
    spectra:
        The spectra of :func:`energy_spectra`, with their variances.
    lines:
        The fits of :func:`fit_elastic_lines`.
    title:
        The title of the figure.
    window:
        The smallest and the largest ``E / Ei`` of the arrival windows, as
        :func:`.arrival_window`.
    shift:
        How far each sub-pulse is lifted above the one below, as a fraction of the highest bin
        of all.

    Returns
    -------
    :
        The figure.
    """
    spectra = spectra.transpose(["subpulse", "energy_transfer"])
    transfer = sc.midpoints(spectra.coords["energy_transfer"].to(unit="meV")).values
    energies = spectra.coords["incident_energy"].to(unit="meV").values
    colors = _subpulse_colors(spectra.sizes["subpulse"])
    step = shift * float(spectra.values.max())
    levels = np.argsort(np.argsort(energies))  # 0 for the lowest Ei, at the bottom
    highest = energies.max()

    with plt.rc_context({"font.size": FONT_SIZE}):  # every text of the figure
        fig, ax = plt.subplots(figsize=(13, 12), layout="constrained")
        fig.suptitle(f"{title}\n{bins_title({'E': spectra.coords['energy_transfer']})}")
        for i, (values, variances, ei, color, level) in enumerate(
            zip(spectra.values, spectra.variances, energies, colors, levels)
        ):
            inside = (transfer >= window[0] * ei) & (transfer <= window[1] * ei)
            x, y, error = transfer[inside], values[inside], np.sqrt(variances[inside])
            height, center, sigma = (
                lines[name].values[i] for name in ("height", "center", "sigma")
            )
            dense = np.linspace(window[0] * ei, window[1] * ei, 4000)
            lift = level * step
            ax.errorbar(  # points, not joined
                x, y + lift, yerr=error, color=color, fmt="o", ms=1.5, elinewidth=0.5, capsize=0
            )
            # E = Ei - Ef: the Gaussian of Ef, centered at Ei - center in E
            ax.plot(dense, gaussian(dense, height, ei - center, sigma) + lift, color=color, lw=1.0)
            fwhm = lines["fwhm"].values[i]
            ax.annotate(
                f"{i}: Ei = {ei:.1f} meV, FWHM = {fwhm:.2f} meV = {100 * fwhm / ei:.1f} % of Ei",
                (window[0] * highest, lift),
                xytext=(4, 5),
                textcoords="offset points",
                va="bottom",
                color=color,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 1},
            )
        ax.set(
            xlim=(window[0] * highest, window[1] * highest),
            xlabel="Energy transfer E = Ei − Ef [meV]",
            ylabel=f"Intensity [{spectra.unit}], each sub-pulse {step:.2g} above the one below",
            title="The elastic line of each sub-pulse, with its Gaussian fit",
        )
        ax.grid(True, alpha=0.3)
    return fig


def q_e_maps(
    events: sc.DataArray,
    pulses: sc.DataArray,
    two_theta: sc.Variable,
    solid_angle: sc.Variable,
    window: tuple[float, float] = DEFAULT_ENERGY_TRANSFER_RATIO,
    energy_step: float = ENERGY_STEP,
    q_step: float = MAP_Q_STEP,
) -> sc.DataArray:
    """Return S(Q, E) of each sub-pulse: its neutrons in |Q| and E, over the coverage of each bin.

    Every sub-pulse has the same bins, so that all the maps span the same range: in E those of
    :func:`energy_spectra`, wide enough for the arrival window of the highest Ei, and in |Q|
    from 0 to the largest |Q| the detector reaches, at the largest 2θ and the highest Ei. In
    each bin, the fraction of the incident neutrons of the sub-pulse scattered into it is
    divided by its coverage from :func:`.q_e_coverage`: the solid angle of every pixel times
    the length in E of its trajectory inside the bin. That is S(Q, E) per steradian and per meV,
    whatever the part of the detector that sees the bin: the normalization of MDNorm, in |Q|.

    Parameters
    ----------
    events:
        The neutrons of :func:`.normalized_events`, in 1/sr, with their ``subpulse``,
        ``energy_transfer``, ``Q`` and the ``solid_angle`` of their pixel.
    pulses:
        The sub-pulses of :func:`.sample_runs.subpulse_timing`.
    two_theta, solid_angle:
        The scattering angle and the solid angle of each pixel the neutrons come from.
    window:
        The smallest and the largest ``E / Ei`` of the arrival windows, as
        :func:`.arrival_window`.
    energy_step, q_step:
        The bin widths, in meV and in 1/Å.

    Returns
    -------
    :
        S(Q, E) per steradian and per meV, along ``subpulse``, ``energy_transfer`` and ``Q``,
        with the bin edges ``energy_transfer`` and ``Q``, in meV and 1/Å, and the
        ``incident_energy`` of each sub-pulse. A bin no trajectory of its sub-pulse crosses,
        which it cannot reach, is NaN.
    """
    energies = pulses.coords["incident_energy"].to(unit="meV")
    highest = float(energies.max().value)
    limit = np.ceil(max(map(abs, window)) * highest)
    energy = bin_edges("energy_transfer", -limit, limit, energy_step, unit="meV")
    # the largest |Q| of all: at the largest 2θ, and the largest kf, of the highest Ei
    largest = float(two_theta.to(unit="deg").max().value)
    reach = float(momentum_transfer(highest, window[0] * highest, largest))
    q = bin_edges("Q", 0.0, float(np.ceil(reach)), q_step, unit="1/Å")
    fractions = events.copy(deep=False)
    fractions.data = events.data * events.coords["solid_angle"]  # of the incident neutrons
    fractions.coords["energy_transfer"] = events.coords["energy_transfer"].to(unit="meV")
    fractions.coords["Q"] = events.coords["Q"].to(unit="1/Å")
    grouped = fractions.group(sc.arange("subpulse", pulses.sizes["subpulse"]))
    scattered = grouped.hist(energy_transfer=energy, Q=q).transpose(
        ["subpulse", "energy_transfer", "Q"]
    )
    coverage = sc.concat(
        [
            q_e_coverage(two_theta, solid_angle, sc.scalar(ei, unit="meV"), energy, q, window)
            for ei in energies.values  # one sub-pulse at a time
        ],
        "subpulse",
    ).to(unit="sr*meV")
    covered = coverage.values > 0
    values = np.full(covered.shape, np.nan)
    values[covered] = scattered.values[covered] / coverage.values[covered]
    per_coverage = sc.scalar(1.0, unit=scattered.unit) / sc.scalar(1.0, unit=coverage.unit)
    maps = sc.DataArray(
        sc.array(dims=scattered.dims, values=values) * per_coverage,
        coords={"energy_transfer": energy, "Q": q},
    )
    maps.coords["incident_energy"] = pulses.coords["incident_energy"].copy()
    return maps


def plot_q_e_maps(maps: sc.DataArray, title: str, columns: int = MAP_COLUMNS) -> Figure:
    """Draw the S(Q,E) map of every sub-pulse, |Q| across and E up, from the highest Ei.

    Every map spans the same range of |Q| and E, that of the bins of :func:`q_e_maps`, and they
    share one linear color scale, turbo, from 0 to the ``MAP_PERCENTILE`` percentile of the bins
    with neutrons, of all the maps. A bin the sub-pulse
    reaches has the color of its value, the bottom of the scale if it is empty; a bin it cannot
    reach, NaN in :func:`q_e_maps`, is white. So the colored area of each map is its coverage:
    its arrival window between the smallest and the largest 2θ of the pixels.

    Parameters
    ----------
    maps:
        The maps of :func:`q_e_maps`.
    title:
        The title of the figure.
    columns:
        The number of maps in a row.

    Returns
    -------
    :
        The figure.
    """
    maps = maps.transpose(["subpulse", "energy_transfer", "Q"])
    q = maps.coords["Q"].to(unit="1/Å").values
    energy = maps.coords["energy_transfer"].to(unit="meV").values
    energies = maps.coords["incident_energy"].to(unit="meV").values
    order = np.argsort(-energies)  # the highest Ei first
    rows = -(-order.size // columns)
    top = float(np.nanpercentile(maps.values[maps.values > 0], MAP_PERCENTILE))
    cmap = plt.get_cmap("turbo").with_extremes(bad="white")
    with plt.rc_context({"font.size": FONT_SIZE}):
        fig, axes = plt.subplots(
            rows,
            columns,
            figsize=(5.5 * columns, 4.5 * rows),
            sharex=True,
            sharey=True,
            layout="constrained",
            squeeze=False,
        )
        fig.suptitle(
            f"{title}\n{bins_title({'|Q|': maps.coords['Q'], 'E': maps.coords['energy_transfer']})}"
        )
        mesh = None
        for ax, i in zip(axes.flat, order):
            ei = float(energies[i])
            mesh = ax.pcolormesh(
                q,
                energy,
                np.ma.masked_invalid(maps.values[i]),  # white where the sub-pulse cannot reach
                cmap=cmap,
                vmin=0.0,
                vmax=top,
            )
            ax.set(
                xlim=(q[0], q[-1]), ylim=(energy[0], energy[-1]), title=f"{i}: Ei = {ei:.1f} meV"
            )
        for ax in axes.flat[order.size :]:
            ax.set_visible(False)
        for ax in axes[-1]:
            ax.set_xlabel("|Q| [1/Å]")
        for ax in axes[:, 0]:
            ax.set_ylabel("E = Ei − Ef [meV]")
        if mesh is not None:
            fig.colorbar(
                mesh,
                ax=axes,
                shrink=0.6,
                extend="max",
                label="S(Q, E) [1/(sr·meV)], normalized by the coverage of each bin",
            )
    return fig


class Reduction(NamedTuple):
    """The reduction of a run of :func:`simulate`, from :func:`reduce_runs`."""

    pulses: sc.DataArray  # the sub-pulses of :func:`.sample_runs.subpulse_timing`
    spectra: sc.DataArray  # of :func:`energy_spectra`
    maps: sc.DataArray  # of :func:`q_e_maps`


def run_names(mode: Mode) -> tuple[str, str]:
    """Return the run folders of a mode: the part to ``sample_pos``, and the part from it.

    Examples
    --------
    >>> run_names("RR")[0]
    'trex_elastic_resolution_incoherent_rr_run_to_sample'
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
        :func:`energy_spectra` and the maps of :func:`q_e_maps`, from the pixels of
        :func:`scattered_pixels`.
    """
    workspace = Workspace.find()
    to_name, from_name = run_names(mode)
    to_sample = Path(to_sample or workspace.run_folder(to_name))
    from_sample = Path(from_sample or workspace.run_folder(from_name))
    # as simulate runs it
    pulses = runs.subpulse_timing_of_run(to_sample, define_instrument(mode))
    binned = scattered_pixels(runs.read_detector(from_sample))
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
        The fits of :func:`fit_elastic_lines`, by mode, with the same sub-pulses.
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


def plot_resolutions(
    lines: dict[Mode, sc.Dataset],
    pulses: dict[Mode, sc.DataArray],
    title: str,
) -> Figure:
    """Draw the FWHM of the elastic line and the flux of each sub-pulse against Ei, in every mode.

    Four panels on one Ei axis, each with one linear y axis: the FWHM, the FWHM as a fraction of Ei,
    the flux at the pre-sample monitor, and the FWHM and the flux of each mode after the first over
    those of the first, the price of the resolution. Each sub-pulse is a point, not joined, in the
    color of its mode, ``MODE_COLORS``; in the last panel an open circle is the FWHM and an open
    square the flux. The FWHM has the error of the fit of :func:`fit_elastic_lines`; the flux has
    none: :func:`.find_patches` sums the intensity of the monitor without its errors, which are
    small from 1e9 rays. The n-th sub-pulse of one mode is the n-th of another, as in
    :func:`resolution_table`.

    Parameters
    ----------
    lines:
        The fits of :func:`fit_elastic_lines`, by mode, with the same sub-pulses.
    pulses:
        The sub-pulses of :func:`.sample_runs.subpulse_timing`, by mode, with their flux.
    title:
        The title of the figure.

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
        absolute.legend(loc="upper left")
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
    ``trex_elastic_resolution_incoherent_m.png``, from :func:`plot_spectra`,
    ``trex_elastic_resolution_incoherent_m_lines.txt``, the table of :func:`lines_table`, and
    ``trex_elastic_resolution_incoherent_m_sqe.png``, from :func:`plot_q_e_maps`; with more than
    one mode, ``trex_elastic_resolution_incoherent_modes.png``, from :func:`plot_resolutions`,
    and ``trex_elastic_resolution_incoherent_modes.txt``, the table of :func:`resolution_table`.

    Parameters
    ----------
    out:
        Folder of the files. If ``None``, ``trex_elastic_resolution_incoherent`` in the results
        folder of :meth:`.Workspace.find`.
    modes:
        The modes to reduce, whose runs :func:`simulate` made.

    Returns
    -------
    :
        The paths of the files.
    """
    folder = cli.output_folder(out, NAME)
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
        )
        files.append(cli.save_figure(fig, folder / f"{NAME}_modes.png"))
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
