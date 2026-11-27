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
"""

import matplotlib.pyplot as plt
import numpy as np
import scipp as sc
from matplotlib.figure import Figure

from mcdgs.diagnostics.diffraction import DIRECT_BEAM_RADIUS
from mcdgs.reduction.binning import bin_edges, bins_title
from mcdgs.reduction.inelastic import DEFAULT_ENERGY_TRANSFER_RATIO

PATTERN_Q_STEP = 0.01  # 1/Å, the bin width in |Q| of the patterns
HEADROOM = 1.05  # the intensity axis reaches this much above the highest bin past the direct beam


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
            unit=scattered.unit / sc.Unit("sr"),
        ),
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
