"""The brilliance transfer of the T-Rex guides, in a window of position and divergence.

The brilliance of a beam is its intensity per area, per solid angle and per wavelength, in
n/s/cm²/sr/Å. Along the beam it can only fall (Liouville's theorem), so the brilliance in a
window over that of the moderator, the brilliance transfer, says how much of what the moderator
could give there the guides deliver: 1 at most.

A :class:`Window` is an area across the beam, centered on its axis, and a divergence of
± its half divergence from the axis, both horizontally and vertically. :func:`main` reads
windows at two points of T-Rex, in the white beam with its choppers parked and the bender out,
as :mod:`.guides_diagnostics` runs it: before the beam monitor, at the exit of the extraction,
and before ``sample_pos``.

The brilliance of the moderator comes from :class:`ModeratorView`: the source alone, and a plane
``REFERENCE_DISTANCE`` in front of it on its axis. A ray that crosses the plane at x within
``REFERENCE_HALF_DIVERGENCE`` of the axis comes from the moderator at about the same x, so the
brilliance there is that of the moderator at x. The thermal face, at x < 0, and the cold face,
at x > 0, are each flat over ``THERMAL_FACE`` and ``COLD_FACE``; the reference is, at each
wavelength, the brighter of the two, averaged over that span.

The figure is on the Guides page of T-Rex, next to those of :mod:`.guides_diagnostics`. Run the
module for it::

    pixi run python -m mcdgs.diagnostics.trex.brilliance_transfer --out docs/instruments/trex/figures/guides
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass, replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from numpy.typing import NDArray

from mcdgs.core import Instrument
from mcdgs.core.components import Arm, ESS_butterfly, Progress_bar
from mcdgs.diagnostics import cli
from mcdgs.instruments.trex.model import Params, TRex
from mcdgs.instruments.trex.sources import source
from mcdgs.mcstas import McStasScriptDiagnostics

NAME = "trex_brilliance_transfer"  # the run folder of T-Rex, and the folder of the figure
REFERENCE_NAME = "trex_moderator_brilliance"  # the run folder of the moderator

# the variables of a ray: position in m, divergence in deg, wavelength in Å and weight in n/s
RAY_VARIABLES = ("x", "y", "dx", "dy", "l", "p")

# the plane of ModeratorView, in front of the moderator, and the rectangle the source aims at on
# it: the width of both faces, and the height of the reference area
REFERENCE_DISTANCE = 1.0  # m
REFERENCE_WIDTH = 0.14  # m
REFERENCE_HEIGHT = 0.01  # m
# deg: through ±0.5° from 1 m, a point of the plane sees ±8.7 mm of the moderator, which keeps
# the view of each face on its flat part
REFERENCE_HALF_DIVERGENCE = 0.5
# m, x on the plane: where each face is flat, from the profile of its brilliance across x
THERMAL_FACE = (-0.07, -0.03)
COLD_FACE = (0.03, 0.07)

# the windows of T-Rex: an area of 1 × 1 cm on the axis, as a sample, and three divergences
AREA_WIDTH = 0.01  # m
AREA_HEIGHT = 0.01  # m
HALF_DIVERGENCES = (0.25, 0.5, 1.0)  # deg
# the points of T-Rex, by the name the figure gives them
POINT_NAMES = {
    "before beam_monitor": "exit of the extraction",
    "before sample_pos": "sample position",
}

# the wavelength bins: 0.1 Å wide, over the white beam, 0.1 to 10 Å
WAVELENGTH_BINS = np.linspace(0.1, 10.0, 100)
# a bin with fewer rays in a window has no brilliance transfer: too few to be read
MIN_RAYS = 20
# the bands of the table that main prints, in Å
BANDS = (
    (0.5, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 3.0), (3.0, 4.0), (4.0, 6.0), (6.0, 8.0), (8.0, 10.0),
)  # fmt: skip

Rays = Mapping[str, NDArray[np.float64]]  # the RAY_VARIABLES of every ray at one point


@dataclass(frozen=True)
class Window:
    """An area across the beam, and a divergence from its axis, both ways.

    Parameters
    ----------
    half_width, half_height:
        Half the size of the area, in m.
    half_divergence:
        The largest angle from the axis, horizontally and vertically, in deg.
    x_center:
        x of the center of the area, in m.
    """

    half_width: float
    half_height: float
    half_divergence: float
    x_center: float = 0.0

    @property
    def area(self) -> float:
        """The area, in cm²."""
        return 4 * self.half_width * self.half_height * 1e4

    @property
    def solid_angle(self) -> float:
        """The solid angle of the divergence, in sr: a square pyramid, 4 asin(sin² a).

        Examples
        --------
        >>> round(Window(0.005, 0.005, 0.5).solid_angle, 9)
        0.00030461
        """
        return 4 * math.asin(math.sin(math.radians(self.half_divergence)) ** 2)

    def contains(self, rays: Rays) -> NDArray[np.bool_]:
        """Tell which rays are in the window, in its area and its divergence."""
        return (
            (np.abs(rays["x"] - self.x_center) <= self.half_width)
            & (np.abs(rays["y"]) <= self.half_height)
            & (np.abs(rays["dx"]) <= self.half_divergence)
            & (np.abs(rays["dy"]) <= self.half_divergence)
        )


def brilliance(
    rays: Rays, window: Window, bins: NDArray[np.float64]
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Compute the brilliance in a window, in each wavelength bin.

    Parameters
    ----------
    rays:
        The rays at the point, with the ``RAY_VARIABLES``, in the frame of the point.
    window:
        The window.
    bins:
        The edges of the wavelength bins, in Å.

    Returns
    -------
    :
        The brilliance, in n/s/cm²/sr/Å, and the number of rays in the window, in each bin.

    Examples
    --------
    >>> rays = {"x": np.zeros(3), "y": np.zeros(3), "dx": np.array([0.0, 0.1, 2.0]),
    ...         "dy": np.zeros(3), "l": np.array([1.5, 1.5, 1.5]), "p": np.ones(3)}
    >>> window = Window(0.005, 0.005, 0.5)
    >>> b, n = brilliance(rays, window, np.array([1.0, 2.0]))
    >>> n.tolist(), bool(np.isclose(b[0], 2 / window.area / window.solid_angle))
    ([2.0], True)
    """
    inside = window.contains(rays)
    intensity, _ = np.histogram(rays["l"][inside], bins=bins, weights=rays["p"][inside])
    counts, _ = np.histogram(rays["l"][inside], bins=bins)
    return intensity / (window.area * window.solid_angle * np.diff(bins)), counts.astype(float)


def moderator_brilliance(
    rays: Rays, bins: NDArray[np.float64]
) -> dict[str, tuple[NDArray[np.float64], NDArray[np.float64]]]:
    """Compute the brilliance of each face of the moderator, from the rays of ModeratorView.

    Parameters
    ----------
    rays:
        The rays on the plane of :class:`ModeratorView`.
    bins:
        The edges of the wavelength bins, in Å.

    Returns
    -------
    :
        The brilliance and the number of rays, as :func:`brilliance` gives them, of
        ``"thermal"``, over ``THERMAL_FACE``, and ``"cold"``, over ``COLD_FACE``.
    """
    faces = {"thermal": THERMAL_FACE, "cold": COLD_FACE}
    return {
        name: brilliance(
            rays,
            Window(
                half_width=(stop - start) / 2,
                half_height=REFERENCE_HEIGHT / 2,
                half_divergence=REFERENCE_HALF_DIVERGENCE,
                x_center=(start + stop) / 2,
            ),
            bins,
        )
        for name, (start, stop) in faces.items()
    }


def reference_brilliance(
    faces: Mapping[str, tuple[NDArray[np.float64], NDArray[np.float64]]],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Take the brighter face of the moderator in each wavelength bin.

    Parameters
    ----------
    faces:
        The brilliance and the number of rays of each face, from :func:`moderator_brilliance`.

    Returns
    -------
    :
        The brilliance of the brighter face, and its number of rays, in each bin.

    Examples
    --------
    >>> faces = {"thermal": (np.array([3.0, 1.0]), np.array([30.0, 10.0])),
    ...          "cold": (np.array([1.0, 2.0]), np.array([20.0, 40.0]))}
    >>> [a.tolist() for a in reference_brilliance(faces)]
    [[3.0, 2.0], [30.0, 40.0]]
    """
    values = np.stack([b for b, _ in faces.values()])
    counts = np.stack([n for _, n in faces.values()])
    brightest = np.argmax(values, axis=0)
    columns = np.arange(values.shape[1])
    return values[brightest, columns], counts[brightest, columns]


def brilliance_transfer(
    value: NDArray[np.float64],
    rays: NDArray[np.float64],
    reference: NDArray[np.float64],
    reference_rays: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Divide a brilliance by that of the moderator, in each bin.

    Parameters
    ----------
    value, rays:
        The brilliance in a window, and its number of rays, from :func:`brilliance`.
    reference, reference_rays:
        The brilliance of the moderator, and its number of rays, from
        :func:`reference_brilliance`.

    Returns
    -------
    :
        The ratio; NaN where either has fewer than ``MIN_RAYS`` rays.

    Examples
    --------
    >>> brilliance_transfer(np.array([1.0, 1.0]), np.array([50.0, 5.0]),
    ...                     np.array([4.0, 4.0]), np.array([100.0, 100.0]))
    array([0.25,  nan])
    """
    enough = (rays >= MIN_RAYS) & (reference_rays >= MIN_RAYS) & (reference > 0)
    return np.divide(value, reference, out=np.full_like(value, np.nan), where=enough)


def band_transfer(
    value: NDArray[np.float64],
    reference: NDArray[np.float64],
    bins: NDArray[np.float64],
    band: tuple[float, float],
) -> float:
    """Divide the brilliance in a band of wavelengths by that of the moderator.

    Parameters
    ----------
    value, reference:
        The brilliance in a window and that of the moderator, in each bin.
    bins:
        The edges of the bins, in Å.
    band:
        The band, in Å: the bins whose centers are in it.

    Returns
    -------
    :
        The ratio of their sums over the band; NaN if the moderator has nothing there.
    """
    centers = (bins[:-1] + bins[1:]) / 2
    in_band = (centers > band[0]) & (centers < band[1])
    total = float(reference[in_band].sum())
    return float(value[in_band].sum()) / total if total > 0 else math.nan


# --8<-- [start:instrument]
@dataclass(frozen=True)
class ModeratorViewParams:
    """Control parameters of :class:`ModeratorView`: there are none."""


class ModeratorView(Instrument):
    """The source of T-Rex alone, and a plane ``REFERENCE_DISTANCE`` in front of it.

    The source aims at a rectangle on the plane, ``REFERENCE_WIDTH`` by ``REFERENCE_HEIGHT``,
    centered on its axis, so every ray crosses the plane within it.

    Parameters
    ----------
    params:
        Control parameters. If ``None``, the defaults of :class:`ModeratorViewParams`.
    """

    params_class = ModeratorViewParams
    name = "trex_moderator_view"  # the McStas instrument name

    origin: Progress_bar
    source: ESS_butterfly
    plane: Arm

    def __init__(self, params: ModeratorViewParams | None = None) -> None:
        super().__init__(params)
        self.origin = Progress_bar(at=[0, 0, 0])
        self.source = replace(
            source(), dist=REFERENCE_DISTANCE, focus_xw=REFERENCE_WIDTH, focus_yh=REFERENCE_HEIGHT
        )
        self.plane = Arm(at=[0, 0, REFERENCE_DISTANCE], relative="source")


def define_moderator_view() -> ModeratorView:
    """Return the instrument of the reference: :class:`ModeratorView`."""
    params = ModeratorViewParams()  # the default control parameters: there are none
    return ModeratorView(params)


def define_instrument() -> TRex:
    """Return the instrument of this diagnostic: T-Rex in the white beam, its choppers parked."""
    # the white beam, 0.1 to 10 Å, every chopper parked on its wide slit (FF), the bender out;
    # the rest as default, as in guides_diagnostics
    params = Params(lambda0=None, rrm=0)
    return TRex(params=params)


# --8<-- [end:instrument]


def plot_brilliance_transfer(
    faces: Mapping[str, tuple[NDArray[np.float64], NDArray[np.float64]]],
    points: Mapping[str, Rays],
    bins: NDArray[np.float64],
    title: str,
) -> Figure:
    """Draw the brilliance, and the brilliance transfer at each point.

    Parameters
    ----------
    faces:
        The brilliance and the number of rays of each face, from :func:`moderator_brilliance`.
    points:
        The rays at each point of ``POINT_NAMES``.
    bins:
        The edges of the wavelength bins, in Å.
    title:
        The title of the figure, such as the running conditions of the instrument.

    Returns
    -------
    :
        Two panels on the same wavelength axis: the brilliance of each face and at each point
        within ±0.5°, on a log scale, and the brilliance transfer of each window of
        ``HALF_DIVERGENCES``, solid at the sample position and dashed at the exit of the
        extraction.
    """
    reference, reference_rays = reference_brilliance(faces)
    fig, (left, right) = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)

    for name, (value, _) in faces.items():
        left.stairs(value, bins, label=f"moderator, {name} face", linewidth=1.5)
    for label, rays in points.items():
        value, _ = brilliance(rays, Window(AREA_WIDTH / 2, AREA_HEIGHT / 2, 0.5), bins)
        left.stairs(value, bins, label=f"{POINT_NAMES[label]}, ±0.5°", linewidth=1.5)
    left.set(ylabel="Brilliance [n/s/cm²/sr/Å]", yscale="log", title="Brilliance")
    # four decades below the brighter face: an empty bin would stretch the axis down to its zero
    left.set_ylim(bottom=float(reference.max()) * 1e-4)

    for color, half_divergence in enumerate(HALF_DIVERGENCES):
        window = Window(AREA_WIDTH / 2, AREA_HEIGHT / 2, half_divergence)
        for label, rays in points.items():
            value, counts = brilliance(rays, window, bins)
            at_sample = label == "before sample_pos"
            right.stairs(
                brilliance_transfer(value, counts, reference, reference_rays),
                bins,
                label=f"{POINT_NAMES[label]}, ±{half_divergence:g}°",
                color=f"C{color}",
                linestyle="-" if at_sample else "--",
                linewidth=1.5,
            )
    right.set(ylabel="Brilliance transfer", title="Over the brighter face of the moderator")
    right.set_ylim(bottom=0)

    for ax in (left, right):
        ax.set(xlabel="Wavelength [Å]")
        ax.set_xlim(bins[0], bins[-1])
        ax.grid(alpha=0.3)
        ax.legend(fontsize="small")
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def _print_bands(
    faces: Mapping[str, tuple[NDArray[np.float64], NDArray[np.float64]]],
    points: Mapping[str, Rays],
    bins: NDArray[np.float64],
) -> None:
    """Print the brilliance transfer of every window at every point, in the ``BANDS``."""
    reference, _ = reference_brilliance(faces)
    print("window | " + " | ".join(f"{start:g}–{stop:g} Å" for start, stop in BANDS))
    for label, rays in points.items():
        for half_divergence in HALF_DIVERGENCES:
            window = Window(AREA_WIDTH / 2, AREA_HEIGHT / 2, half_divergence)
            value, _ = brilliance(rays, window, bins)
            transfers = (band_transfer(value, reference, bins, band) for band in BANDS)
            print(
                f"{POINT_NAMES[label]}, ±{half_divergence:g}° | "
                + " | ".join(f"{t:.3f}" for t in transfers)
            )


def main(
    out: Path | None = None,
    ncount: float = 1e7,
    reference_ncount: float = 4e6,
    seed: int = 1,
) -> Path:
    """Run the moderator alone and T-Rex, and save the brilliance and its transfer.

    The figure is ``trex_brilliance_transfer.png``. The runs go to the run folders
    ``trex_moderator_brilliance`` and ``trex_brilliance_transfer``, replaced on every run. It
    also prints the brilliance transfer in the ``BANDS``.

    Parameters
    ----------
    out:
        Folder of the figure. If ``None``, ``trex_brilliance_transfer`` in the results folder
        of :meth:`.Workspace.find`.
    ncount:
        Number of McStas neutrons of T-Rex, which keeps every ray at both points.
    reference_ncount:
        Number of McStas neutrons of the moderator alone, which keeps every ray on its plane.
    seed:
        Seed of the random numbers.

    Returns
    -------
    :
        The path of the figure.
    """
    # --8<-- [start:reference]
    moderator = McStasScriptDiagnostics(define_moderator_view())
    moderator.add_point(before="plane")
    moderator.add_view("x", "l", bins=100)  # every run needs a view
    moderator.run(REFERENCE_NAME, ncount=reference_ncount, seed=seed)
    plane = moderator.rays(*RAY_VARIABLES)["before plane"]
    faces = moderator_brilliance(plane, WAVELENGTH_BINS)
    # --8<-- [end:reference]

    # --8<-- [start:trex]
    diag = McStasScriptDiagnostics(define_instrument())
    for label in POINT_NAMES:  # before each component: in its frame, on its axis
        diag.add_point(before=label.removeprefix("before "))
    diag.add_view("x", "y", bins=100)
    diag.run(NAME, ncount=ncount, seed=seed)
    points = diag.rays(*RAY_VARIABLES)
    # --8<-- [end:trex]

    _print_bands(faces, points, WAVELENGTH_BINS)
    title = diag.simulation.instrument.running_conditions()
    fig = plot_brilliance_transfer(faces, points, WAVELENGTH_BINS, title)
    return cli.save_figure(fig, cli.output_folder(out, NAME) / f"{NAME}.png", dpi=100)


if __name__ == "__main__":
    main(out=cli.parse_out(__doc__))
