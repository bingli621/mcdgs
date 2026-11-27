"""Diffraction from a powder, reduced to I(|Q|) and drawn, for any instrument.

The powder-diffraction diagnostic of each instrument runs McStas, finds its sub-pulses, reads its
detector and converts every neutron in the arrival window of its sub-pulse to its momentum
transfer with :func:`.normalized_events`, as :mod:`.diagnostics.diffraction` does for a crystal.
A powder scatters the same way about the beam, so its pattern depends on |Q| alone. This module
does the rest:

- :func:`powder_patterns` reduces the neutrons of each sub-pulse to its intensity in the elastic
  |Q| of their pixel, per incident neutron and per steradian, over the solid angle of the pixels
  that see each bin.
- :func:`plot_powder_patterns` draws the patterns of all the sub-pulses over each other, on a
  linear scale, with their errors.
- :func:`detector_views` reduces the neutrons of each sub-pulse to its intensity on the detector,
  in the horizontal angle θ and the height of their pixel, per incident neutron and per
  steradian, over the solid angle of the pixels in each bin.
- :func:`plot_detector_views` draws the view of each sub-pulse on its own panel, on one log
  color scale.
"""

import matplotlib.pyplot as plt
import numpy as np
import scipp as sc
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure

from mcdgs.diagnostics.diffraction import DIRECT_BEAM_RADIUS
from mcdgs.reduction.binning import bin_edges, bins_title
from mcdgs.reduction.inelastic import DEFAULT_ENERGY_TRANSFER_RATIO

PATTERN_Q_STEP = 0.01  # 1/Å, the bin width in |Q| of the patterns
HEADROOM = 1.05  # the intensity axis reaches this much above the highest bin past the direct beam
# the bins of the detector views, at least the pixels of both detectors, so that no bin inside a
# detector is empty: the banana detector of T-Rex has 0.625° by 16 mm pixels, the panels of
# CSPEC tubes 0.42° apart, 10 mm pixels
THETA_STEP = 1.0  # deg
HEIGHT_STEP = 0.02  # m
DIRECT_BEAM_THETA = 5.0  # deg: the color scale of the detector views leaves out |θ| below this
# the color scale of the detector views spans this many decades below the strongest line, down to
# the incoherent scattering of the powder, about 1e-5 of it
VIEW_LOG_DECADES = 6
DETECTOR_COLUMNS = 2  # the detector views of the sub-pulses are drawn this many to a row
DETECTOR_PANEL_WIDTH = 9.0  # inches, of each detector view


def powder_patterns(
    events: sc.DataArray,
    pulses: sc.DataArray,
    two_theta: sc.Variable,
    solid_angle: sc.Variable,
    step: float = PATTERN_Q_STEP,
) -> sc.DataArray:
    """Return the intensity in the elastic |Q| of each sub-pulse, over the coverage of each bin.

    Each neutron of :func:`.normalized_events` times the solid angle of its pixel is the
    fraction of the incident neutrons of its sub-pulse it stands for, with the variance of its
    square: a bin has the error √Σw². It is put at the elastic |Q| of its pixel, ``2 ki sin θ``
    with the ki of its sub-pulse and the 2θ of the center of the pixel, not at the |Q| of its
    time of flight: in each bin, the sum is divided by its coverage, the solid angle of the
    pixels whose elastic |Q| falls in it, so a neutron and its pixel fall in the same bin. That
    is the intensity per incident neutron and per steradian, whatever the part of the detector
    that sees the bin. With the |Q| of the time of flight, a neutron could fall in a bin between
    two columns of pixels, which none of them covers. The powder scatters elastically, so the two
    differ only by the resolution.

    Parameters
    ----------
    events:
        The neutrons of :func:`.normalized_events`, in 1/sr, with their ``subpulse``, and the
        ``solid_angle`` and ``two_theta`` of their pixel.
    pulses:
        The sub-pulses, with their ``incident_wavelength`` and ``incident_energy``.
    two_theta, solid_angle:
        The scattering angle and the solid angle of every pixel of the detector.
    step:
        The bin width in |Q|, in 1/Å.

    Returns
    -------
    :
        The intensity per incident neutron, at the monitor of :func:`.normalized_events`, and
        per steradian, with its variances, along
        ``subpulse`` and ``Q``, with the bin edges ``Q`` in 1/Å, from 0 past the largest elastic
        |Q| of the pixels, and the ``incident_wavelength`` and ``incident_energy`` of each
        sub-pulse. A bin no pixel of its sub-pulse sees is NaN.
    """
    ki = 2 * np.pi / pulses.coords["incident_wavelength"].to(unit="Å").values
    index = events.coords["subpulse"]
    fractions = events.copy(deep=False)
    fractions.data = events.data * events.coords["solid_angle"]
    fractions.variances = fractions.values**2  # a neutron of weight w counts with variance w²
    # the elastic |Q| of the pixel of each neutron, 2 ki sin θ, with the ki of its sub-pulse
    half_angle = np.radians(events.coords["two_theta"].to(unit="deg").values) / 2
    fractions.coords["Q"] = sc.array(
        dims=events.dims, values=2 * ki[index.values] * np.sin(half_angle), unit="1/Å"
    )
    # that of every pixel, for every sub-pulse
    pixel_half_angle = np.radians(two_theta.to(unit="deg").values) / 2
    pixel_q = 2 * ki[:, np.newaxis] * np.sin(pixel_half_angle)[np.newaxis, :]
    top = float(np.ceil(pixel_q.max() * 10) / 10)
    edges = bin_edges("Q", 0.0, top, step, unit="1/Å")
    n = pulses.sizes["subpulse"]
    groups = sc.arange("subpulse", n, unit=index.unit, dtype=index.dtype)
    scattered = fractions.group(groups).hist(Q=edges).transpose(["subpulse", "Q"])
    weights = solid_angle.to(unit="sr").values
    coverage = np.stack([np.histogram(q, edges.values, weights=weights)[0] for q in pixel_q])
    covered = coverage > 0
    values = np.full(coverage.shape, np.nan)
    variances = np.full(coverage.shape, np.nan)
    values[covered] = scattered.values[covered] / coverage[covered]
    variances[covered] = scattered.variances[covered] / coverage[covered] ** 2
    patterns = sc.DataArray(
        sc.array(
            dims=["subpulse", "Q"],
            values=values,
            variances=variances,
            unit=scattered.unit,
        )
        / sc.scalar(1.0, unit="sr"),
        coords={"Q": edges},
    )
    for name in ("incident_wavelength", "incident_energy"):
        patterns.coords[name] = pulses.coords[name].copy()
    return patterns


def plot_powder_patterns(
    patterns: sc.DataArray,
    title: str,
    monitor: str,
    window: tuple[float, float] = DEFAULT_ENERGY_TRANSFER_RATIO,
    direct_beam_radius: float = DIRECT_BEAM_RADIUS,
) -> Figure:
    """Draw the intensity in |Q| of every sub-pulse over each other, with its errors.

    The scale is linear in intensity, from 0 to just above the highest bin past the direct beam,
    ``|Q| > direct_beam_radius``: the direct beam that passes the beamstop, far above the Bragg
    peaks, runs off the top. Each sub-pulse has its own color, of turbo, from the shortest
    wavelength, and is drawn where its pixels see the bin, each bin as a point with its error,
    joined to the next. The bins are written under the title, as :func:`.bins_title`.

    Parameters
    ----------
    patterns:
        The patterns of :func:`powder_patterns`.
    title:
        The title of the figure.
    monitor:
        The monitor whose intensity of each sub-pulse the patterns are normalized by, as
        :func:`.normalized_events` does, such as ``"pre-sample monitor"``: it names the unit of
        the intensity axis, per neutron of the sub-pulse there.
    window:
        The smallest and the largest ``E / Ei`` of the arrival windows, as
        :func:`.arrival_window`, which the patterns sum over.
    direct_beam_radius:
        The |Q| of the direct beam, in 1/Å, that the intensity axis leaves out.

    Returns
    -------
    :
        The figure.
    """
    q = sc.midpoints(patterns.coords["Q"]).values
    wavelengths = patterns.coords["incident_wavelength"].to(unit="Å").values
    energies = patterns.coords["incident_energy"].to(unit="meV").values
    cmap = plt.get_cmap("turbo")
    colors = [cmap(x) for x in np.linspace(0.05, 0.95, patterns.sizes["subpulse"])]

    fig, ax = plt.subplots(figsize=(14, 8), layout="constrained")
    bins = bins_title({"|Q|": patterns.coords["Q"], "E/Ei": bin_edges("E/Ei", *window)})
    fig.suptitle(f"{title}\n{bins}")
    for i, (values, variances, color) in enumerate(
        zip(patterns.values, patterns.variances, colors)
    ):
        seen = np.isfinite(values)  # where the pixels of the sub-pulse see the bin
        ax.errorbar(
            q[seen],
            values[seen],
            yerr=np.sqrt(variances[seen]),
            color=color,
            fmt="o-",  # points, joined
            ms=1.5,
            lw=0.8,
            elinewidth=0.5,
            capsize=0,
            label=f"{i}: λi = {wavelengths[i]:.3f} Å, Ei = {energies[i]:.2f} meV",
        )
    past_beam = patterns.values[:, q > direct_beam_radius]
    ax.set(
        ylim=(0, HEADROOM * float(np.nanmax(past_beam))),
        xlim=(0, q[-1]),
        xlabel="|Q| [1/Å]",
        ylabel=f"Intensity [{patterns.unit} per neutron at the {monitor}]",
    )
    ax.grid(True, alpha=0.3)
    ax.legend(title="sub-pulse", fontsize="small", loc="upper right")  # above the high-|Q| peaks
    return fig


def horizontal_angle(position: sc.Variable, sample_position: sc.Variable) -> sc.Variable:
    """Return the horizontal angle θ of each position about the vertical through the sample.

    θ is 0 along the beam, +z, and positive toward +x, as the ``theta`` of the banana detector
    of T-Rex and the :func:`.panel_angle` of the panels of CSPEC.

    Parameters
    ----------
    position:
        The positions, such as of the pixels.
    sample_position:
        The position of the sample.

    Returns
    -------
    :
        The angle in deg, from -180 to 180, along the dimensions of ``position``.
    """
    offset = position - sample_position
    return sc.atan2(y=offset.fields.x, x=offset.fields.z).to(unit="deg")


def covering_range(values: np.ndarray, step: float) -> tuple[float, float]:
    """Return the first and the last edge of bins of ``step`` aligned on 0, covering ``values``.

    Parameters
    ----------
    values:
        The values to cover.
    step:
        The bin width.

    Returns
    -------
    :
        The edges, multiples of ``step``, the last above the largest value.
    """
    return float(np.floor(values.min() / step) * step), float(
        (np.floor(values.max() / step) + 1) * step
    )


def detector_views(
    events: sc.DataArray,
    pulses: sc.DataArray,
    position: sc.Variable,
    sample_position: sc.Variable,
    solid_angle: sc.Variable,
    theta_step: float = THETA_STEP,
    height_step: float = HEIGHT_STEP,
) -> sc.DataArray:
    """Return the intensity of each sub-pulse on the detector, in θ and height.

    Each neutron is put at the horizontal angle θ of :func:`horizontal_angle` and the height of
    the center of its pixel, as a detector shows it, unrolled about the vertical through the
    sample. As in :func:`powder_patterns`, each neutron of :func:`.normalized_events` times the
    solid angle of its pixel is the fraction of the incident neutrons of its sub-pulse it stands
    for, with the variance of its square, and in each bin the sum is divided by its coverage, the
    solid angle of the pixels whose center falls in it: the intensity per incident neutron and per
    steradian, whatever the number of pixels in the bin. A bin wider than a pixel has at least
    one, inside the detector. The radius of the detector is the median horizontal distance of the
    pixels from the sample, to draw the views to the scale of the detector.

    Parameters
    ----------
    events:
        The neutrons of :func:`.normalized_events`, in 1/sr, with their ``subpulse``, the
        ``solid_angle`` and ``position`` of their pixel.
    pulses:
        The sub-pulses, with their ``incident_wavelength`` and ``incident_energy``.
    position, solid_angle:
        The position and the solid angle of every pixel of the detector.
    sample_position:
        The position of the sample.
    theta_step:
        The bin width in θ, in deg.
    height_step:
        The bin width in height, in m.

    Returns
    -------
    :
        The intensity per incident neutron, at the monitor of :func:`.normalized_events`, and
        per steradian, with its variances, along ``subpulse``, ``height`` and ``theta``, with the
        bin edges ``theta`` in deg and ``height`` in m, covering the pixels, and the
        ``incident_wavelength`` and ``incident_energy`` of each sub-pulse, and the ``radius`` of
        the detector, in m. A bin without a pixel is NaN.
    """
    pixel_theta = horizontal_angle(position, sample_position).values
    pixel_offset = position - sample_position
    pixel_height = pixel_offset.fields.y.to(unit="m").values
    radius = np.median(np.hypot(pixel_offset.fields.x.values, pixel_offset.fields.z.values))
    edges = {
        "theta": bin_edges("theta", *covering_range(pixel_theta, theta_step), theta_step, "deg"),
        "height": bin_edges("height", *covering_range(pixel_height, height_step), height_step, "m"),
    }
    fractions = events.copy(deep=False)
    fractions.data = events.data * events.coords["solid_angle"]
    fractions.variances = fractions.values**2  # a neutron of weight w counts with variance w²
    offset = events.coords["position"] - sample_position
    fractions.coords["theta"] = horizontal_angle(events.coords["position"], sample_position)
    fractions.coords["height"] = offset.fields.y.to(unit="m")
    index = events.coords["subpulse"]
    groups = sc.arange("subpulse", pulses.sizes["subpulse"], unit=index.unit, dtype=index.dtype)
    scattered = fractions.group(groups).hist(height=edges["height"], theta=edges["theta"])
    coverage, _, _ = np.histogram2d(
        pixel_height,
        pixel_theta,
        bins=[edges["height"].values, edges["theta"].values],
        weights=solid_angle.to(unit="sr").values,
    )
    covered = np.broadcast_to(coverage > 0, scattered.shape)
    values = np.divide(
        scattered.values, coverage, out=np.full(scattered.shape, np.nan), where=covered
    )
    variances = np.divide(
        scattered.variances, coverage**2, out=np.full(scattered.shape, np.nan), where=covered
    )
    views = sc.DataArray(
        sc.array(
            dims=["subpulse", "height", "theta"],
            values=values,
            variances=variances,
            unit=scattered.unit,
        )
        / sc.scalar(1.0, unit="sr"),
        coords=edges,
    )
    for name in ("incident_wavelength", "incident_energy"):
        views.coords[name] = pulses.coords[name].copy()
    views.coords["radius"] = sc.scalar(float(radius), unit=pixel_offset.unit).to(unit="m")
    return views


def plot_detector_views(
    views: sc.DataArray,
    title: str,
    monitor: str,
    window: tuple[float, float] = DEFAULT_ENERGY_TRANSFER_RATIO,
    log_decades: float = VIEW_LOG_DECADES,
    direct_beam_theta: float = DIRECT_BEAM_THETA,
    headroom: float = 1.0,
    columns: int = DETECTOR_COLUMNS,
) -> Figure:
    """Draw the detector view of each sub-pulse on its own panel, θ across and the height up.

    The panels are drawn to the scale of the detector: a degree of θ is as long across as its arc,
    ``radius · π / 180``, is high, so that near the beam, where a Bragg cone crosses the
    cylinder of the detector in a circle of radius ``radius · 2θ``, the ring is round.

    θ runs from high to low, left to right. All the panels share the axes and one log color scale, in turbo, with the bins without a
    neutron or a pixel in white. Its bottom is ``log_decades`` below the highest bin of any
    sub-pulse past the direct beam, ``|θ| > direct_beam_theta``, and its top ``headroom`` times
    that bin: the direct beam past the beamstop takes the color of the top. The bins are written under the title, as
    :func:`.bins_title`.

    Parameters
    ----------
    views:
        The views of :func:`detector_views`.
    title:
        The title of the figure.
    monitor:
        The monitor whose intensity of each sub-pulse the views are normalized by, as
        :func:`plot_powder_patterns`.
    window:
        The smallest and the largest ``E / Ei`` of the arrival windows, as
        :func:`.arrival_window`, which the views sum over.
    log_decades:
        How many decades below its top the color scale starts.
    direct_beam_theta:
        The |θ| of the direct beam, in deg, that the color scale leaves out.
    headroom:
        The top of the color scale, as a multiple of the highest bin past the direct beam. Above
        1, the strongest lines are no longer at the red end.
    columns:
        How many panels to a row.

    Returns
    -------
    :
        The figure.
    """
    theta = views.coords["theta"].values
    height = views.coords["height"].values
    past_beam = np.abs(sc.midpoints(views.coords["theta"]).values) > direct_beam_theta
    highest = float(np.nanmax(views.values[..., past_beam]))
    norm = LogNorm(vmin=highest / 10**log_decades, vmax=highest * headroom, clip=True)
    cmap = plt.get_cmap("turbo").with_extremes(bad="white")
    wavelengths = views.coords["incident_wavelength"].to(unit="Å").values
    energies = views.coords["incident_energy"].to(unit="meV").values
    # a meter of height is as long as a meter of arc, 180 / (π · radius) degrees of θ
    aspect = 180 / (np.pi * views.coords["radius"].to(unit="m").value)

    n = views.sizes["subpulse"]
    rows = -(-n // columns)  # ceiling division
    # panels as high as their aspect, plus room for their titles
    span = aspect * (height[-1] - height[0]) / (theta[-1] - theta[0])
    row_height = DETECTOR_PANEL_WIDTH * span + 0.6
    fig, axes = plt.subplots(
        rows,
        columns,
        figsize=(DETECTOR_PANEL_WIDTH * columns + 1.5, row_height * rows + 1.2),
        sharex=True,
        sharey=True,
        squeeze=False,
        layout="constrained",
    )
    bins = bins_title(
        {
            "θ": views.coords["theta"],
            "height": views.coords["height"],
            "E/Ei": bin_edges("E/Ei", *window),
        }
    )
    fig.suptitle(f"{title}\n{bins}")
    for i, ax in enumerate(axes.flat):
        if i >= n:
            ax.set_visible(False)  # the panels left over in the last row
            continue
        values = views["subpulse", i].values
        ax.pcolormesh(
            theta,
            height,
            np.ma.masked_where(~(values > 0), values),  # no neutron or no pixel: white
            cmap=cmap,
            norm=norm,
        )
        ax.set_aspect(aspect)  # to the scale of the detector: the rings are round
        ax.set_xlim(theta[-1], theta[0])  # θ from high to low
        ax.set_title(
            f"sub-pulse {i}\nλi = {wavelengths[i]:.3f} Å, Ei = {energies[i]:.2f} meV",
            fontsize="small",
        )
        ax.grid(True, alpha=0.3)
    for ax in axes[-1]:
        ax.set_xlabel("θ [deg], from the beam")
    for ax in axes[:, 0]:
        ax.set_ylabel("Height [m]")
    fig.colorbar(
        ScalarMappable(norm=norm, cmap=cmap),  # the scale that all the panels share
        ax=axes,
        shrink=0.6,
        label=f"Intensity [{views.unit} per neutron at the {monitor}]",
        extend="max",
    )
    return fig
