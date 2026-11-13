"""The banana detector of T-Rex lit by vanadium in the white beam.

:class:`TRexWithVanadium` is :class:`.TRex` with the vanadium cylinder of :func:`.vanadium` on
``sample_pos``, and every chopper parked with a slit on the beam, so the sample is in the white
beam: ``lambda0`` is ``None``, so the source simulates its whole band, ``L_MIN`` to ``L_MAX``,
at every time of the pulse. Vanadium scatters it elastically
and evenly in every direction, so each pixel of the banana detector counts in proportion to the
solid angle it spans, less what the sample absorbs on the way out.

:func:`detector_image` sums the neutrons of each pixel, on the grid of the detector: the angle
2θ of the pixel in the horizontal plane, from the beam, and its height from the beam axis.
:func:`plot_detector_image` draws it.

:func:`pixel_intensity` gives the intensity of each pixel against its number, with the solid
angle of the pixel from :func:`.pixel_solid_angle`, and :func:`plot_pixel_intensity` draws
both: the solid angle times one factor for the whole detector, the intensity per steradian of
:func:`intensity_per_solid_angle`, follows the intensity.

:func:`main` runs McStas, from the source to the detector, then draws both figures. Run the
module for them::

    pixi run python -m mcdgs.diagnostics.trex.solid_angles
"""

import math
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipp as sc
from matplotlib.figure import Figure

from mcdgs.core.params import InstrParam
from mcdgs.diagnostics.cli import output_folder, parse_out, save_figure
from mcdgs.diagnostics.trex.elastic_plane import pixel_solid_angle
from mcdgs.instruments.trex.choppers import DISKS, slit_angle
from mcdgs.instruments.trex.detectors import (
    BANANA_CENTER_Y,
    BANANA_HEIGHT,
    BANANA_RADIUS,
    BANANA_THETA_BINS,
    BANANA_THETA_LIMITS,
    BANANA_Y_BINS,
)
from mcdgs.instruments.trex.model import Params, TRex
from mcdgs.instruments.trex.monitors import beamstop_monitor, pre_sample_monitor, sample_monitor
from mcdgs.instruments.trex.sources import wavelength_band
from mcdgs.mcstas.nexus_io import read_events
from mcdgs.mcstas.simulation import Simulation
from mcdgs.samples import vanadium

DETECTOR = "banana_detector"
SAMPLE = "sample_pos"
SOURCE = "pre_sample_monitor"  # on the beam just before the sample: gives its direction there
NAME = "trex_solid_angles"
RUN = f"{NAME}_run"  # the run folder
FOCUS_MARGIN = 1.0  # deg, the focusing window of the sample is this much wider than the detector
SAMPLE_RADIUS = 0.01  # m, a thin cylinder: little of what it scatters is absorbed on the way out
SAMPLE_SPLIT = 10  # McStas SPLIT on the sample: each neutron reaching it is repeated this often
WHITE_BEAM = Params(lambda0=None)  # no central wavelength: the band is not limited
# s, the width of a time bin of the monitors near the sample: the white beam takes 0.4 s to
# pass them, which the bins of T-Rex, made for one sub-pulse, would split into 80000
MONITOR_TIME_RESOLUTION = 1e-3
DIRECT_BEAM_ANGLE = 5.0  # deg: the direct beam, far above the rest, is left out of the scale


def detector_window() -> dict[str, float]:
    """Return the focusing window that covers the banana detector, seen from the sample.

    The window is centered on the middle of the detector in the horizontal plane, and reaches
    as far above that plane as below it, to the lower edge of the detector, which is the further
    one. McStas turns a window centered above or below the plane as a whole, which would leave
    the corners of the detector outside. It is ``FOCUS_MARGIN`` wider on every side: the sample
    is not a point.

    Returns
    -------
    :
        ``target_x``, ``target_z``, ``focus_aw`` and ``focus_ah`` of the sample, in the frame
        of ``sample_pos``.

    Examples
    --------
    >>> window = detector_window()
    >>> round(window["focus_aw"], 2), round(window["focus_ah"], 2)
    (182.0, 51.96)
    """
    theta_start, theta_stop = BANANA_THETA_LIMITS
    # the elevation of the edge of the detector furthest from the horizontal plane
    reach = max(abs(BANANA_CENTER_Y + sign * BANANA_HEIGHT / 2) for sign in (-1, 1))
    elevation = math.degrees(math.atan2(reach, BANANA_RADIUS))
    azimuth = math.radians((theta_start + theta_stop) / 2)
    return {
        "target_x": math.sin(azimuth),
        "target_z": math.cos(azimuth),
        "focus_aw": theta_stop - theta_start + 2 * FOCUS_MARGIN,
        "focus_ah": 2 * (elevation + FOCUS_MARGIN),
    }


class TRexWithVanadium(TRex):
    """T-Rex with a vanadium sample in the white beam: every chopper is parked open.

    The sample is that of :func:`.vanadium`, a cylinder 10 mm high, thinned to a radius of
    ``SAMPLE_RADIUS``, which scatters towards the window of
    :func:`detector_window` only, split ``SAMPLE_SPLIT`` times. Each disk stands still
    (``nu = 0``) with the slit of the mode on the beam, as :func:`.slit_angle`: the wide ones in
    ``"FF"``. There are no instrument parameters ``nu_<name>`` and ``delay_<name>``.

    The time and wavelength monitors near the sample (pre-sample, beamstop and sample) span the
    whole band, with time bins of ``MONITOR_TIME_RESOLUTION``.

    The source is the plain one: the masked source spends its rays where turning choppers let
    neutrons through.

    Parameters
    ----------
    params:
        Control parameters. ``lambda0`` gives the band of the beam, ``None`` for the white
        beam, and ``mode`` the slit each disk is parked on.
    """

    name = "trex_vanadium"  # the McStas instrument name

    def __init__(self, params: Params = WHITE_BEAM) -> None:
        sample = vanadium(
            name="sample", radius=SAMPLE_RADIUS, split=SAMPLE_SPLIT, **detector_window()
        )
        super().__init__(sample=sample, params=params)
        band = wavelength_band(params.lambda0)
        for monitor in (pre_sample_monitor, beamstop_monitor, sample_monitor):
            wide = monitor(lambda_limits=band, time_resolution=MONITOR_TIME_RESOLUTION)
            setattr(self, monitor.__name__, wide)
        for name in DISKS:
            angle = slit_angle(name, params.mode)
            parked = replace(getattr(self, name), nu=0, delay=0, park_angle=angle)
            setattr(self, name, parked)

    def chopper_instr_params(self) -> list[InstrParam]:
        """Return no instrument parameters for the choppers: they are parked."""
        return []


def detector_image(binned: sc.DataArray) -> sc.DataArray:
    """Sum the neutrons of each pixel, on the grid of the banana detector.

    Parameters
    ----------
    binned:
        The neutrons by pixel, as :func:`.read_events`, in n/s, with the pixel ``position``, the
        ``sample_position`` and the ``source_position``, in the lab frame. The beam at the
        sample runs from the source to the sample: the curved guides have turned it from the
        z axis of the lab.

    Returns
    -------
    :
        The intensity of each pixel, in n/s, along ``height``, the height of the pixel from the
        beam axis, in m, and ``two_theta``, its angle in the horizontal plane from the beam, in
        deg, positive to the left of the beam. The bin edges are those of the pixels:
        ``BANANA_Y_BINS`` by ``BANANA_THETA_BINS``.
    """
    relative = binned.coords["position"] - binned.coords["sample_position"]
    pixels = binned.bins.sum()
    pixels.coords["two_theta"] = horizontal_angle(binned)
    pixels.coords["height"] = relative.fields.y.to(unit="m")
    theta_start, theta_stop = BANANA_THETA_LIMITS
    bottom = BANANA_CENTER_Y - BANANA_HEIGHT / 2
    return pixels.hist(
        height=sc.linspace("height", bottom, bottom + BANANA_HEIGHT, BANANA_Y_BINS + 1, unit="m"),
        two_theta=sc.linspace(
            "two_theta", theta_start, theta_stop, BANANA_THETA_BINS + 1, unit="deg", dtype=float
        ),
    )


def plot_detector_image(
    image: sc.DataArray, title: str, direct_beam_angle: float = DIRECT_BEAM_ANGLE
) -> Figure:
    """Draw the intensity of each pixel of the banana detector, unrolled.

    The color scale is turbo, linear from 0 to the highest pixel beyond the direct beam,
    ``|2θ| > direct_beam_angle``: the direct beam that passes the beamstop, a few pixels far
    above all others, would leave the rest at the bottom of the scale, so it takes the color of
    the top instead.

    Parameters
    ----------
    image:
        The intensity of each pixel, as :func:`detector_image`.
    title:
        The title of the figure.
    direct_beam_angle:
        The |2θ| of the direct beam, in deg, that the color scale leaves out.

    Returns
    -------
    :
        The figure.
    """
    two_theta = sc.midpoints(image.coords["two_theta"]).values
    vmax = float(image.values[:, np.abs(two_theta) > direct_beam_angle].max())
    fig, ax = plt.subplots(figsize=(12, 5), layout="constrained")
    fig.suptitle(title)
    mesh = ax.pcolormesh(
        image.coords["two_theta"].values,
        image.coords["height"].values,
        image.values,
        cmap="turbo",
        vmin=0,
        vmax=vmax,
    )
    ax.set(
        xlabel="2θ [deg], in the horizontal plane",
        ylabel="Height [m], from the beam axis",
    )
    fig.colorbar(mesh, ax=ax, label=f"Intensity [{image.unit} per pixel]", extend="max")
    return fig


def horizontal_angle(binned: sc.DataArray) -> sc.Variable:
    """Return the angle of each pixel in the horizontal plane, from the beam at the sample.

    Parameters
    ----------
    binned:
        The neutrons by pixel, as :func:`.read_events`, with the pixel ``position``, the
        ``sample_position`` and the ``source_position``, in the lab frame. The beam at the
        sample runs from the source to the sample: the curved guides have turned it from the
        z axis of the lab.

    Returns
    -------
    :
        The angle of each pixel, in deg, positive to the left of the beam.
    """
    relative = binned.coords["position"] - binned.coords["sample_position"]
    beam = binned.coords["sample_position"] - binned.coords["source_position"]
    return (
        sc.atan2(y=relative.fields.x, x=relative.fields.z)
        - sc.atan2(y=beam.fields.x, x=beam.fields.z)
    ).to(unit="deg")


def pixel_intensity(binned: sc.DataArray) -> sc.DataArray:
    """Sum the neutrons of each pixel, against its number, with its solid angle.

    The variance of a pixel is the sum of the squares of the weights of its neutrons, as McStas
    takes it for a monitor.

    Parameters
    ----------
    binned:
        The neutrons by pixel, as :func:`.read_events`, in n/s.

    Returns
    -------
    :
        The intensity of each pixel that was hit, in n/s, with its variance, along
        ``pixel_id``, in the order of the pixel numbers, with the coordinates ``solid_angle``, in sr, from
        :func:`.pixel_solid_angle`, and ``two_theta``, from :func:`horizontal_angle`.
    """
    pixels = binned.bins.sum()
    pixels.variances = (binned.bins.data * binned.bins.data).bins.sum().values
    pixels.coords["solid_angle"] = pixel_solid_angle(binned)
    pixels.coords["two_theta"] = horizontal_angle(binned)
    return sc.sort(pixels, "pixel_id")


def intensity_per_solid_angle(
    pixels: sc.DataArray, direct_beam_angle: float = DIRECT_BEAM_ANGLE
) -> sc.Variable:
    """Return the one factor between the intensity of the pixels and their solid angle.

    It is the intensity of all the pixels beyond the direct beam, ``|2θ| > direct_beam_angle``,
    over the solid angle they span together.

    Parameters
    ----------
    pixels:
        The intensity of each pixel, as :func:`pixel_intensity`.
    direct_beam_angle:
        The |2θ| of the direct beam, in deg, left out.

    Returns
    -------
    :
        The intensity per solid angle, in n/s/sr, with its variance.
    """
    outside = np.abs(pixels.coords["two_theta"].values) > direct_beam_angle
    (dim,) = pixels.dims
    kept = pixels[dim, np.flatnonzero(outside)]
    return kept.data.sum() / kept.coords["solid_angle"].sum()


def plot_pixel_intensity(
    pixels: sc.DataArray, title: str, direct_beam_angle: float = DIRECT_BEAM_ANGLE
) -> Figure:
    """Draw the intensity of each pixel against its number, with its solid angle over it.

    Each pixel has an error bar of one standard deviation, from its variance. The solid angle is times the factor of :func:`intensity_per_solid_angle`. The intensity
    axis ends at the highest pixel beyond the direct beam, as in :func:`plot_detector_image`.

    Parameters
    ----------
    pixels:
        The intensity of each pixel, as :func:`pixel_intensity`.
    title:
        The title of the figure.
    direct_beam_angle:
        The |2θ| of the direct beam, in deg, left out of the factor and of the intensity axis.

    Returns
    -------
    :
        The figure.
    """
    factor = intensity_per_solid_angle(pixels, direct_beam_angle)
    pixel_id = pixels.coords["pixel_id"].values
    outside = np.abs(pixels.coords["two_theta"].values) > direct_beam_angle

    fig, ax = plt.subplots(figsize=(12, 5), layout="constrained")
    fig.suptitle(title)
    ax.errorbar(
        pixel_id,
        pixels.values,
        yerr=np.sqrt(pixels.variances),
        fmt=".",
        markersize=1,
        elinewidth=0.2,
        color="tab:blue",
        label="Simulated, ± 1 standard deviation",
    )
    ax.plot(
        pixel_id,
        (sc.values(factor) * pixels.coords["solid_angle"]).values,
        ".",
        markersize=1,
        color="tab:red",
        zorder=3,  # over the error bars
        label=f"Solid angle × {factor.value:.4g} n/s/sr",
    )
    ax.set(
        xlabel="Pixel number",
        ylabel=f"Intensity [{pixels.unit} per pixel]",
        xlim=(pixel_id.min(), pixel_id.max()),
        ylim=(0, 1.05 * pixels.values[outside].max()),
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right", markerscale=10)
    return fig


def main(out: Path | None = None, ncount: float = 1e8, seed: int = 1) -> list[Path]:
    """Run :class:`TRexWithVanadium` in McStas, and save the figures of the detector.

    The run goes from the source to the detector, into ``RUN`` in :func:`.default_output_dir`,
    and replaces the folder of an earlier one. The figures are
    ``trex_solid_angles_detector.png``, from :func:`plot_detector_image`, and
    ``trex_solid_angles_pixels.png``, from :func:`plot_pixel_intensity`.

    Parameters
    ----------
    out:
        Folder of the figures. If ``None``, the folder of the run, next to its ``mccode.h5``.
    ncount:
        Number of McStas rays.
    seed:
        Seed of the random numbers.

    Returns
    -------
    :
        The paths of the figures, in that order.
    """
    instrument = TRexWithVanadium()  # the white beam
    run = Simulation(instrument).run(RUN, ncount=ncount, seed=seed).parent

    binned = read_events(run, DETECTOR, SOURCE, SAMPLE)
    image = detector_image(binned)
    hit = np.count_nonzero(image.values)
    print(
        f"  {binned.bins.size().sum().value:.0f} neutrons on {hit} of {image.values.size} pixels,"
        f" {image.sum().value:.4g} {image.unit}"
    )
    l_min, l_max = instrument.source.Lmin, instrument.source.Lmax
    conditions = (
        f"{instrument.name} · {instrument.params.mode} · choppers parked open · "
        f"λ = {l_min:.3f} to {l_max:.3f} Å"
    )
    folder = output_folder(out, run)
    title = "T-Rex, vanadium in the white beam: the banana detector"
    image_figure = plot_detector_image(image, f"{title}\n{conditions}")
    title = "T-Rex, vanadium in the white beam: the pixels of the banana detector"
    pixels_figure = plot_pixel_intensity(pixel_intensity(binned), f"{title}\n{conditions}")
    return [
        save_figure(image_figure, folder / f"{NAME}_detector.png"),
        save_figure(pixels_figure, folder / f"{NAME}_pixels.png"),
    ]


if __name__ == "__main__":
    main(out=parse_out(__doc__))
