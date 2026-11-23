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

:func:`main` runs it in two parts, with :class:`.Simulation`, as
:mod:`.diagnostics.trex.elastic_plane` does: from the source to ``sample_pos``, which dumps the
beam there into a beam file, then from that beam through the vanadium to the detector.

The reduction is that of :mod:`.diagnostics.trex.elastic_plane`:

1. :func:`.elastic_plane.subpulse_timing` finds the sub-pulses at the pre-sample and beamstop
   monitors of the first part, with their incident energy and the time they reach the sample.
2. :func:`.read_events` reads the neutrons of the banana detector of the second part, by pixel,
   and :func:`scattered_pixels` leaves out the pixels within ``MIN_TWO_THETA`` of the beam: the
   direct beam that passes the beamstop.
3. :func:`.elastic_plane.normalized_events` keeps the neutrons inside the arrival window of a
   sub-pulse, gives them their energy transfer, and divides each by the solid angle of its pixel
   and by the flux of its sub-pulse.
4. :func:`energy_spectra` sums them into a histogram of the energy transfer of each sub-pulse,
   averaged over the solid angle of the pixels: per incident neutron, per steradian and per meV.

:func:`fit_elastic_lines` fits a Gaussian to the elastic line of each sub-pulse, against the
final energy Ef, and gives its FWHM, ``FWHM_PER_SIGMA`` times its standard deviation σ.
:func:`plot_spectra` draws the spectra with their errors and their Gaussians, against E.
:func:`q_e_maps` histograms the neutrons of each sub-pulse in |Q| and E as well, and
:func:`plot_q_e_maps` draws the maps of all of them in one figure.

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
from mcdgs.diagnostics.trex import elastic_plane as plane
from mcdgs.diagnostics.trex.solid_angles import detector_window
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.mcstas import Simulation, Workspace, read_events
from mcdgs.reduction.binning import bin_edges
from mcdgs.reduction.inelastic import DEFAULT_ENERGY_TRANSFER_RATIO
from mcdgs.samples import vanadium

NAME = "trex_elastic_resolution_incoherent"
TO_SAMPLE = f"{NAME}_run_to_sample"  # the run folder of the first part, source -> sample_pos
FROM_SAMPLE = f"{NAME}_run_from_sample"  # that of the second part, sample_pos -> detector
MIN_TWO_THETA = 5.0  # deg: the pixels closer to the beam see the direct beam past the beamstop
ENERGY_STEP = 0.05  # meV, the bin width of the spectra
FWHM_PER_SIGMA = 2 * np.sqrt(2 * np.log(2))  # the full width at half maximum of a Gaussian
# the Gaussian is fitted to the bins of at least this fraction of the highest one: further out,
# a bin holds few neutrons, and the variance of their weights is no error for it
FIT_FRACTION = 0.2
SHIFT = 0.4  # each sub-pulse is lifted above the one below by this much of the highest bin
FONT_SIZE = 14  # pt, of every text of the figure
MAP_Q_STEP = 0.05  # 1/Å, the bin width in |Q| of the S(Q,E) maps; in E, ENERGY_STEP
# the unit of the maps, as scipp cannot print meV and 1/Å together but in SI units
MAP_UNIT = "1/(sr·meV·Å⁻¹)"
MAP_COLUMNS = 3  # the maps of the sub-pulses, in rows of this many
# meV Å², ħ²/2m of the neutron: E = HBAR2_OVER_2M k²
HBAR2_OVER_2M = float(
    (sc.constants.hbar**2 / (2 * sc.constants.neutron_mass)).to(unit="meV*Å**2").value
)


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


def define_instrument() -> TRexElasticResolution:
    """Return the instrument of this diagnostic: :class:`TRexElasticResolution` at its defaults."""
    params = Params()  # the default control parameters
    return TRexElasticResolution(params=params)


# --8<-- [end:instrument]


def two_theta(binned: sc.DataArray) -> sc.Variable:
    """Return the scattering angle of each pixel, from the beam along ``BEAM_DIRECTION``.

    Parameters
    ----------
    binned:
        The neutrons by pixel, as :func:`.read_events`, with the pixel ``position`` and the
        ``sample_position``, in the beam frame of ``sample_pos``.

    Returns
    -------
    :
        The angle between the beam and the direction from the sample to each pixel, in deg,
        along the pixel dimension.
    """
    scattered = binned.coords["position"] - binned.coords["sample_position"]
    beam = sc.vector(plane.BEAM_DIRECTION)
    cosine = sc.dot(scattered, beam) / sc.norm(scattered)
    return sc.acos(sc.values(cosine)).to(unit="deg")


def scattered_pixels(binned: sc.DataArray, min_two_theta: float = MIN_TWO_THETA) -> sc.DataArray:
    """Leave out the pixels close to the beam, which see the direct beam past the beamstop.

    Parameters
    ----------
    binned:
        The neutrons by pixel, as :func:`.read_events`.
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

    Each neutron of :func:`.elastic_plane.normalized_events` times the solid angle of its pixel
    is the fraction of the incident neutrons of its sub-pulse it stands for. Summed over the
    pixels and divided by their solid angle, it is the intensity averaged over the detector.
    Each neutron of weight w counts with the variance w², so a bin has the error √Σw².

    Parameters
    ----------
    events:
        The neutrons of :func:`.elastic_plane.normalized_events`, in 1/sr, with their
        ``subpulse``, ``energy_transfer`` and the ``solid_angle`` of their pixel.
    pulses:
        The sub-pulses of :func:`.elastic_plane.subpulse_timing`.
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
        fig.suptitle(title)
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


def momentum_transfer(
    incident_energy: float, energy_transfer: np.ndarray, two_theta: float
) -> np.ndarray:
    """Return |Q| of a neutron of an incident energy and energy transfer, scattered by 2θ.

    ``Q² = ki² + kf² - 2 ki kf cos 2θ``, with ``E = HBAR2_OVER_2M k²``.

    Parameters
    ----------
    incident_energy:
        Ei, in meV.
    energy_transfer:
        E = Ei - Ef, in meV, below Ei.
    two_theta:
        The scattering angle, in deg.

    Returns
    -------
    :
        |Q| for each energy transfer, in 1/Å.

    Examples
    --------
    Elastic and scattered straight back, Q is 2 ki:

    >>> ki = np.sqrt(20.0 / HBAR2_OVER_2M)
    >>> bool(np.isclose(momentum_transfer(20.0, np.array([0.0]), 180.0)[0], 2 * ki))
    True
    """
    ki = np.sqrt(incident_energy / HBAR2_OVER_2M)
    kf = np.sqrt((incident_energy - energy_transfer) / HBAR2_OVER_2M)
    return np.sqrt(ki**2 + kf**2 - 2 * ki * kf * np.cos(np.radians(two_theta)))


def q_e_maps(
    events: sc.DataArray,
    pulses: sc.DataArray,
    solid_angle: sc.Variable,
    two_theta_limits: tuple[float, float],
    window: tuple[float, float] = DEFAULT_ENERGY_TRANSFER_RATIO,
    energy_step: float = ENERGY_STEP,
    q_step: float = MAP_Q_STEP,
) -> sc.DataArray:
    """Histogram the neutrons of each sub-pulse in |Q| and E, as :func:`energy_spectra` does in E.

    Every sub-pulse has the same bins, so that all the maps span the same range: in E those of
    :func:`energy_spectra`, wide enough for the arrival window of the highest Ei, and in |Q|
    from 0 to the largest |Q| the detector reaches, at the largest 2θ and the highest Ei. The
    intensity is that of :func:`energy_spectra` per bin of |Q| too, so a map summed over |Q|
    gives the spectrum of its sub-pulse. It is averaged over the solid angle of the whole
    detector, not normalized by the coverage of each bin.

    Parameters
    ----------
    events:
        The neutrons of :func:`.elastic_plane.normalized_events`, in 1/sr, with their
        ``subpulse``, ``energy_transfer``, ``Q`` and the ``solid_angle`` of their pixel.
    pulses:
        The sub-pulses of :func:`.elastic_plane.subpulse_timing`.
    solid_angle:
        The solid angle of all the pixels the neutrons come from, in sr.
    two_theta_limits:
        The smallest and the largest 2θ of those pixels, in deg.
    window:
        The smallest and the largest ``E / Ei`` of the arrival windows, as
        :func:`.arrival_window`.
    energy_step, q_step:
        The bin widths, in meV and in 1/Å.

    Returns
    -------
    :
        The intensity per incident neutron, per steradian, per meV and per 1/Å
        (``MAP_UNIT``), along ``subpulse``, ``energy_transfer`` and ``Q``, with the bin edges
        ``energy_transfer`` and ``Q``, in meV and 1/Å, and the ``incident_energy`` of each
        sub-pulse.
    """
    energies = pulses.coords["incident_energy"].to(unit="meV").values
    highest = float(energies.max())
    limit = np.ceil(max(map(abs, window)) * highest)
    energy = bin_edges("energy_transfer", -limit, limit, energy_step, unit="meV")
    # the largest |Q| of all: at the largest 2θ, and the largest kf, of the highest Ei
    reach = momentum_transfer(highest, np.array([window[0] * highest]), two_theta_limits[1])[0]
    q = bin_edges("Q", 0.0, float(np.ceil(reach)), q_step, unit="1/Å")
    fractions = events.copy(deep=False)
    fractions.data = events.data * events.coords["solid_angle"]  # of the incident neutrons
    fractions.coords["energy_transfer"] = events.coords["energy_transfer"].to(unit="meV")
    fractions.coords["Q"] = events.coords["Q"].to(unit="1/Å")
    grouped = fractions.group(sc.arange("subpulse", pulses.sizes["subpulse"]))
    maps = grouped.hist(energy_transfer=energy, Q=q) / solid_angle.to(unit="sr")
    maps = maps / sc.scalar(energy_step, unit="meV") / sc.scalar(q_step, unit="1/Å")
    maps.coords["incident_energy"] = pulses.coords["incident_energy"].copy()
    return maps


def plot_q_e_maps(
    maps: sc.DataArray,
    two_theta_limits: tuple[float, float],
    title: str,
    window: tuple[float, float] = DEFAULT_ENERGY_TRANSFER_RATIO,
    columns: int = MAP_COLUMNS,
) -> Figure:
    """Draw the S(Q,E) map of every sub-pulse, |Q| across and E up, from the highest Ei.

    Every map spans the same range of |Q| and E, that of the bins of :func:`q_e_maps`, and they
    share one linear color scale, turbo, from 0 to the highest bin of all, with the empty bins
    in white. The dashed lines are the |Q| of the smallest and the largest 2θ of the detector at
    each E of the arrival window of the sub-pulse, ``window`` times its Ei: the detector reaches
    the bins between them only.

    Parameters
    ----------
    maps:
        The maps of :func:`q_e_maps`.
    two_theta_limits:
        The smallest and the largest 2θ of the pixels, in deg.
    title:
        The title of the figure.
    window:
        The smallest and the largest ``E / Ei`` of the arrival windows, as
        :func:`.arrival_window`.
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
    top = float(maps.values.max())
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
        fig.suptitle(title)
        mesh = None
        for ax, i in zip(axes.flat, order):
            ei = float(energies[i])
            mesh = ax.pcolormesh(
                q,
                energy,
                np.ma.masked_less_equal(maps.values[i], 0),  # empty bins in white
                cmap=cmap,
                vmin=0.0,
                vmax=top,
            )
            own = np.linspace(window[0] * ei, window[1] * ei, 400)
            for two_theta in two_theta_limits:
                ax.plot(momentum_transfer(ei, own, two_theta), own, "k--", lw=0.8)
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
                label=f"S(Q, E) [{MAP_UNIT}], not normalized by the coverage of each bin",
            )
    return fig


class Reduction(NamedTuple):
    """The reduction of a run of :func:`simulate`, from :func:`reduce_runs`."""

    pulses: sc.DataArray  # the sub-pulses of :func:`.elastic_plane.subpulse_timing`
    spectra: sc.DataArray  # of :func:`energy_spectra`
    maps: sc.DataArray  # of :func:`q_e_maps`
    two_theta_limits: tuple[float, float]  # deg, the smallest and the largest 2θ of the pixels


def simulate(ncount: float = 1e9, seed: int = 1) -> tuple[Path, Path]:
    """Run :class:`TRexElasticResolution` in McStas, in two parts split at ``sample_pos``.

    The first part runs from the source to ``sample_pos`` into ``TO_SAMPLE``, with the monitors
    of the incident energy and the beam file; the second from that beam through the vanadium to
    the detector, into ``FROM_SAMPLE``, both in the runs folder of :meth:`.Workspace.find`. Each
    replaces the folder of an earlier run.

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
    sim = Simulation(define_instrument())
    to_sample = sim.run_to_sample(folder_name=TO_SAMPLE, ncount=ncount, seed=seed)
    # sample_pos -> detector, from the beam file, in the frame of sample_pos: the beam frame
    from_sample = sim.run_from_sample(folder_name=FROM_SAMPLE, seed=seed + 1)
    return to_sample.folder, from_sample.folder


def reduce_runs(
    to_sample: str | Path | None = None, from_sample: str | Path | None = None
) -> Reduction:
    """Reduce the two parts of a McStas run of :func:`simulate`, without running McStas.

    Parameters
    ----------
    to_sample, from_sample:
        The folders of the first and the second part. If ``None``, ``TO_SAMPLE`` and
        ``FROM_SAMPLE`` in the runs folder of :meth:`.Workspace.find`, where :func:`simulate`
        writes.

    Returns
    -------
    :
        The sub-pulses of :func:`.elastic_plane.subpulse_timing`, the spectra of
        :func:`energy_spectra` and the maps of :func:`q_e_maps`, from the pixels of
        :func:`scattered_pixels`, and the smallest and the largest 2θ of those pixels.
    """
    workspace = Workspace.find()
    to_sample = Path(to_sample or workspace.run_folder(TO_SAMPLE))
    from_sample = Path(from_sample or workspace.run_folder(FROM_SAMPLE))
    pulses = plane.subpulse_timing_of_run(to_sample)
    binned = read_events(
        from_sample, plane.DETECTOR, "source", plane.SAMPLE, beam_direction=plane.BEAM_DIRECTION
    )
    binned = scattered_pixels(binned)
    events = plane.normalized_events(binned, pulses)
    solid_angle = plane.pixel_solid_angle(binned).sum()
    angles = two_theta(binned)
    limits = (float(angles.min().value), float(angles.max().value))
    return Reduction(
        pulses,
        energy_spectra(events, pulses, solid_angle),
        q_e_maps(events, pulses, solid_angle, limits),
        limits,
    )


def analyze(
    out: Path | None = None,
    to_sample: str | Path | None = None,
    from_sample: str | Path | None = None,
) -> list[Path]:
    """Reduce the two parts of a McStas run, and save the spectra and the fits of their lines.

    It runs no McStas, so it can run again on the same runs; :func:`main` runs them first. The
    files are ``trex_elastic_resolution_incoherent.png``, from :func:`plot_spectra`,
    ``trex_elastic_resolution_incoherent_lines.txt``, the table of :func:`lines_table`, and
    ``trex_elastic_resolution_incoherent_sqe.png``, from :func:`plot_q_e_maps`.

    Parameters
    ----------
    out:
        Folder of the files. If ``None``, ``trex_elastic_resolution_incoherent`` in the results
        folder of :meth:`.Workspace.find`.
    to_sample, from_sample:
        The folders of the two parts, as :func:`reduce_runs`.

    Returns
    -------
    :
        The paths of the files.
    """
    pulses, spectra, maps, limits = reduce_runs(to_sample, from_sample)
    folder = cli.output_folder(out, NAME)
    conditions = define_instrument().running_conditions()
    lines = fit_elastic_lines(spectra)
    table = lines_table(lines)
    print(table, end="")
    files = [cli.save_text(table, folder / f"{NAME}_lines.txt")]
    fig = plot_spectra(
        spectra,
        lines,
        f"T-Rex, vanadium: the energy transfer of {pulses.sizes['subpulse']} sub-pulses\n"
        f"{conditions}",
    )
    files.append(cli.save_figure(fig, folder / f"{NAME}.png"))
    fig = plot_q_e_maps(maps, limits, f"T-Rex, vanadium: S(Q, E) of each sub-pulse\n{conditions}")
    files.append(cli.save_figure(fig, folder / f"{NAME}_sqe.png"))
    return files


def main(out: Path | None = None, ncount: float = 1e9, seed: int = 1) -> list[Path]:
    """Run McStas with :func:`simulate`, then reduce and draw it with :func:`analyze`.

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
    to_sample, from_sample = simulate(ncount=ncount, seed=seed)
    return analyze(out, to_sample, from_sample)


if __name__ == "__main__":
    analyze(out=cli.parse_out(__doc__))  # McStas runs only with main, a few minutes
