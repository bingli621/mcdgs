"""The elastic resolution from vanadium, reduced and drawn, for any instrument.

Vanadium scatters incoherently and elastically, so every neutron it scatters has an energy
transfer of zero: the width of its elastic line is the energy resolution of the instrument. The
elastic resolution diagnostic of each instrument runs McStas, finds its sub-pulses, reads its
detector and gives each neutron its energy transfer with :func:`.normalized_events`; this module
does the rest:

- :func:`two_theta` gives the scattering angle of each pixel, and :func:`scattered_pixels` leaves
  out the pixels close to the beam.
- :func:`energy_spectra` sums the neutrons into a histogram of the energy transfer of each
  sub-pulse, averaged over the solid angle of the pixels: per incident neutron, per steradian and
  per meV. Like :func:`q_e_maps`, it divides each neutron by its kf/ki, with
  :func:`scattered_fractions`, so that a line has the shape of S(E), not of (kf/ki) S(E).
- :func:`fit_elastic_lines` fits a Gaussian to the elastic line of each sub-pulse, against the
  final energy Ef, and gives its FWHM, ``FWHM_PER_SIGMA`` times its standard deviation σ;
  :func:`lines_table` lists them, and :func:`plot_spectra` draws the spectra with their errors and
  their Gaussians, against E.
- :func:`q_e_maps` histograms the neutrons of each sub-pulse in |Q| and E as well, divided by
  their kf/ki, into S(Q, E), and :func:`plot_q_e_maps` draws the maps of all of them in one
  figure.

The sub-pulses are a data array along ``subpulse``, with their flux and the coordinates
``incident_energy``, ``incident_wavelength`` and ``time_at_sample``, as
:func:`.sample_runs.subpulse_timing` gives them for T-Rex.
"""

import matplotlib.pyplot as plt
import numpy as np
import scipp as sc
from matplotlib.figure import Figure
from scipy.optimize import curve_fit

from mcdgs.diagnostics import cli
from mcdgs.reduction.binning import bin_edges, bins_title
from mcdgs.reduction.inelastic import DEFAULT_ENERGY_TRANSFER_RATIO, kf_over_ki
from mcdgs.reduction.mdnorm import momentum_transfer, q_e_coverage

# the incident beam in the frame of sample_pos, the frame of a run from the sample: its z axis
BEAM_DIRECTION = (0.0, 0.0, 1.0)
ENERGY_STEP = 0.05  # meV, the bin width of the spectra
FWHM_PER_SIGMA = 2 * np.sqrt(2 * np.log(2))  # the full width at half maximum of a Gaussian
# the Gaussian is fitted to the bins of at least this fraction of the highest one: further out,
# a bin holds few neutrons, and the variance of their weights is no error for it
FIT_FRACTION = 0.2
SHIFT = 0.4  # each sub-pulse is lifted above the one below by this much of the highest bin
FONT_SIZE = 14  # pt, of every text of the figure
MAP_Q_STEP = 0.05  # 1/Å, the bin width in |Q| of the S(Q,E) maps; in E, ENERGY_STEP
MAP_COLUMNS = 5  # the maps of the sub-pulses, in rows of this many
# the top of the color scale of the maps, as a percentile of the bins with neutrons: a bin on the
# edge of the coverage, with little of it, can stand far above the rest
MAP_PERCENTILE = 99.9


def two_theta(binned: sc.DataArray) -> sc.Variable:
    """Return the scattering angle of each pixel, from the beam along ``BEAM_DIRECTION``.

    Parameters
    ----------
    binned:
        The neutrons by pixel, as :func:`.read_events`, with the pixel ``position`` and the
        ``sample_position``, in the frame of ``sample_pos``.

    Returns
    -------
    :
        The angle between the beam and the direction from the sample to each pixel, in deg,
        along the pixel dimension.
    """
    scattered = binned.coords["position"] - binned.coords["sample_position"]
    beam = sc.vector(BEAM_DIRECTION)
    cosine = sc.dot(scattered, beam) / sc.norm(scattered)
    return sc.acos(sc.values(cosine)).to(unit="deg")


def scattered_pixels(binned: sc.DataArray, min_two_theta: float) -> sc.DataArray:
    """Leave out the pixels close to the beam, such as those that see the direct beam.

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


def scattered_fractions(events: sc.DataArray, pulses: sc.DataArray) -> sc.DataArray:
    """Return the fraction of the incident neutrons each neutron stands for, divided by kf/ki.

    Each neutron of :func:`.normalized_events` times the solid angle of its pixel is the fraction
    of the incident neutrons of its sub-pulse it stands for. The cross section holds kf/ki on top
    of S(Q, E): above 1 where the neutron gains energy, below where it loses it. Each neutron is
    divided by its kf/ki, :func:`.kf_over_ki`, with the Ei of its sub-pulse, so that what is
    summed is S(Q, E).

    Parameters
    ----------
    events:
        The neutrons of :func:`.normalized_events`, in 1/sr, with their ``subpulse``,
        ``energy_transfer`` and the ``solid_angle`` of their pixel.
    pulses:
        The sub-pulses, with their ``incident_energy``.

    Returns
    -------
    :
        The neutrons, dimensionless, with their ``energy_transfer`` in meV.
    """
    fractions = events.copy(deep=False)
    fractions.coords["energy_transfer"] = events.coords["energy_transfer"].to(unit="meV")
    incident = pulses.coords["incident_energy"].to(unit="meV").values
    ratio = kf_over_ki(
        sc.array(dims=events.dims, values=incident[events.coords["subpulse"].values], unit="meV"),
        fractions.coords["energy_transfer"],
    )
    fractions.data = events.data * events.coords["solid_angle"] / ratio
    return fractions


def energy_spectra(
    events: sc.DataArray,
    pulses: sc.DataArray,
    solid_angle: sc.Variable,
    step: float = ENERGY_STEP,
) -> sc.DataArray:
    """Sum the neutrons of each sub-pulse into a histogram of their energy transfer.

    Each neutron of :func:`.normalized_events` times the solid angle of its pixel
    is the fraction of the incident neutrons of its sub-pulse it stands for, divided by its kf/ki
    with :func:`scattered_fractions`. Summed over the pixels and divided by their solid angle, it
    is S(E) averaged over the detector. Each neutron of weight w counts with the variance w², so a
    bin has the error √Σw².

    Parameters
    ----------
    events:
        The neutrons of :func:`.normalized_events`, in 1/sr, with their
        ``subpulse``, ``energy_transfer`` and the ``solid_angle`` of their pixel.
    pulses:
        The sub-pulses, with their ``incident_energy`` and ``incident_wavelength``.
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
    fractions = scattered_fractions(events, pulses)
    fractions.variances = fractions.values**2  # a neutron of weight w counts with variance w²
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
    from 0 to the largest |Q| the detector reaches, at the largest 2θ and the highest Ei. Each
    neutron is divided by its kf/ki, :func:`.kf_over_ki`, which the cross section holds on top of
    S(Q, E): it is 1 at E = 0, above 1 where the neutron gains energy and below where it loses
    it. In each bin, the fraction of the incident neutrons of the sub-pulse scattered into it is
    then divided by its coverage from :func:`.q_e_coverage`: the solid angle of every pixel times
    the length in E of its trajectory inside the bin. That is S(Q, E) per steradian and per meV,
    whatever the part of the detector that sees the bin: the normalization of MDNorm, in |Q|.

    Parameters
    ----------
    events:
        The neutrons of :func:`.normalized_events`, in 1/sr, with their ``subpulse``,
        ``energy_transfer``, ``Q`` and the ``solid_angle`` of their pixel.
    pulses:
        The sub-pulses, with their ``incident_energy`` and ``incident_wavelength``.
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
    fractions = scattered_fractions(events, pulses)
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


def plot_q_e_maps(
    maps: sc.DataArray,
    title: str,
    columns: int = MAP_COLUMNS,
    percentile: float = MAP_PERCENTILE,
) -> Figure:
    """Draw the S(Q,E) map of every sub-pulse, |Q| across and E up, from the highest Ei.

    Every map spans the same range of |Q| and E, that of the bins of :func:`q_e_maps`, and they
    share one linear color scale, turbo, from 0 to the ``percentile`` percentile of the bins with
    neutrons, of all the maps. A bin the sub-pulse reaches has the color of its value, the bottom of the scale if it is empty; a bin it cannot
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
    percentile:
        The top of the color scale, as a percentile of the bins with neutrons. A lower one
        saturates more bins, and shows the weaker ones in more of the scale.

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
    top = float(np.nanpercentile(maps.values[maps.values > 0], percentile))
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
