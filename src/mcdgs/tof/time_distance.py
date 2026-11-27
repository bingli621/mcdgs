"""The time-distance diagram of any :class:`.Instrument`, from its tof model.

:func:`time_distance_diagram` builds the tof model of the instrument (:func:`.build_tof`), runs
it and draws tof's time-distance diagram: each neutron is a ray from the moderator, through the
openings of the choppers, to the detectors, colored by its wavelength, and the neutrons a
chopper stops are drawn in gray up to it. Unlike the acceptance diagram (:mod:`.acceptance`), it
follows neutrons, so it shows the shape of the pulse, with Monte Carlo noise.

With ``between``, the diagram zooms in on the flight from one component of the tof model to
another, such as from a monitor before the sample to the detector, where an inelastic sample
(:mod:`.sample`) spreads the arrival times: in distance from the one to the other, and in time
from the earliest neutron at the first to the latest at the second, of those that reach the
second.
"""

import numpy as np
import tof
from matplotlib.figure import Figure

from mcdgs.core import Instrument
from mcdgs.tof.builder import build_tof

ZOOM_MARGIN = 0.05  # of the span of a zoomed axis, added on either side


def time_distance_diagram(
    instrument: Instrument,
    neutrons: int = 1_000_000,
    *,
    seed: int | None = None,
    n_pulses: int = 1,
    visible_rays: int = 2000,
    blocked_rays: int = 5000,
    figsize: tuple[float, float] = (10, 6),
    cmap: str = "turbo",
    between: tuple[str, str] | None = None,
    title: str | None = None,
    color_by: str | None = None,
    color_headroom: float = 0.0,
    time_percentile: float = 100.0,
) -> Figure:
    """Build the tof model of an instrument, run it and draw its time-distance diagram.

    The title is the instrument's :meth:`.Instrument.running_conditions`, under ``title`` if
    one is given. tof draws the time in µs.

    Parameters
    ----------
    instrument:
        The instrument, at its current parameters.
    neutrons:
        Number of tof neutrons per pulse.
    seed:
        Seed of the random numbers, of the source and of the choice of rays. If ``None``, every
        run differs.
    n_pulses:
        Number of pulses, one source period apart.
    visible_rays:
        Most rays drawn, per pulse, of the neutrons that reach the last detector.
    blocked_rays:
        Most rays drawn in gray, per pulse, of the neutrons a chopper stops. 0 draws none.
    figsize:
        Size of the figure, in inches.
    cmap:
        Colormap of the wavelength.
    between:
        The names of two components of the tof model, the nearer first, to zoom in on the flight
        between them (:func:`zoom_between`). If ``None``, the whole instrument.
    title:
        A line above the running conditions.
    color_by:
        The name of a component of the tof model whose wavelengths span the color scale, from
        the shortest to the longest of the neutrons that reach the furthest component. tof
        colors each stretch of a ray by its wavelength where it starts, so the stretches before
        an inelastic sample take their incident wavelength. With the component before the
        sample, each sub-pulse gets its own color there, and the stretches after the sample past
        either end of the scale take its end colors. If ``None``, the wavelengths at the
        furthest component, as tof draws it.
    color_headroom:
        With ``color_by``, Å added to the top of the color scale, so that more of the
        wavelengths after an inelastic sample, lengthened by its energy loss, keep a color of
        their own.
    time_percentile:
        With ``between``, where the time axis ends, as :func:`zoom_between`.

    Returns
    -------
    :
        The figure.

    Raises
    ------
    ValueError
        If the instrument has no ``ESS_butterfly`` or ``Polygon_ESS_butterfly``, or more than
        one; if a chopper is parked; or as :func:`zoom_between`.
    """
    model = build_tof(instrument, neutrons=neutrons, seed=seed, n_pulses=n_pulses)
    result = model.run()
    vmin = vmax = None
    if color_by is not None:
        vmin, vmax = wavelength_range(result, color_by, furthest_component(model))
        vmax += color_headroom
    conditions = instrument.running_conditions()
    plot = result.plot(
        visible_rays=visible_rays,
        blocked_rays=blocked_rays,
        figsize=figsize,
        cmap=cmap,
        seed=seed,
        vmin=vmin,
        vmax=vmax,
        title=conditions if title is None else f"{title}\n{conditions}",
    )
    if between is not None:
        plot.ax.set(**zoom_between(result, *between, percentile=time_percentile))
        for text in plot.ax.texts:  # the names of the components and the pulses outside it
            text.set_clip_on(True)
    plot.fig.tight_layout()
    return plot.fig


def furthest_component(model: tof.Model) -> str:  # type: ignore
    """Return the name of the component of a tof model furthest from the source."""
    return max(model.components.values(), key=lambda c: c.distance.value).name


def wavelength_range(result: tof.Result, name: str, reaching: str) -> tuple[float, float]:  # type: ignore
    """Return the shortest and the longest wavelength at a component of a tof run.

    Parameters
    ----------
    result:
        The tof run.
    name:
        The component whose wavelengths are taken.
    reaching:
        The component the neutrons must reach, such as the furthest: those it records unmasked.

    Returns
    -------
    :
        The wavelengths, in Å, as tof draws them.

    Raises
    ------
    ValueError
        If a name is not a component of the run, or no neutron reaches ``reaching``.
    """
    for component in (name, reaching):
        try:
            result[component]
        except KeyError:
            raise ValueError(f"{component!r} is not a component of the tof model") from None
    through = ~result[reaching].data.masks["blocked_by_others"].values
    wavelengths = result[name].data.coords["wavelength"].to(unit="Å").values[through]
    if not np.isfinite(wavelengths).any():
        raise ValueError(f"no neutron reaches {reaching!r}")
    return float(np.nanmin(wavelengths)), float(np.nanmax(wavelengths))


def zoom_between(
    result: tof.Result,  # type: ignore
    first: str,
    second: str,
    percentile: float = 100.0,
) -> dict[str, tuple[float, float]]:
    """Return the limits of a time-distance diagram on the flight between two components.

    In distance, from the first to the second; in time, from the earliest arrival at the first
    to the latest at the second, of the neutrons that reach the second, over every pulse, or to
    the ``percentile`` percentile of their arrivals at the second. Each axis has ``ZOOM_MARGIN``
    of its span on either side.

    Parameters
    ----------
    result:
        The tof run.
    first, second:
        The names of its components, the nearer first.
    percentile:
        Where the time axis ends, as a percentile of the arrivals at the second: below 100, the
        few slowest neutrons, such as those an inelastic sample leaves with almost no energy,
        run past it.

    Returns
    -------
    :
        ``ylim``, the distances in m, and ``xlim``, the times in µs, as tof draws them.

    Raises
    ------
    ValueError
        If a name is not a component of the run, if the first is not nearer the source than the
        second, or if no neutron reaches the second.
    """
    readings = {}
    for name in (first, second):
        try:
            readings[name] = result[name]
        except KeyError:
            raise ValueError(f"{name!r} is not a component of the tof model") from None
    near, far = (readings[name].distance.to(unit="m").value for name in (first, second))
    if near >= far:
        raise ValueError(
            f"{first!r}, {near:g} m, is not nearer the source than {second!r}, {far:g} m"
        )
    # the neutrons that reach the second, and their times at the first and the second
    reached = readings[second].data
    through = ~reached.masks["blocked_by_others"].values
    start = readings[first].data.coords["toa"].to(unit="µs").values[through]
    end = reached.coords["toa"].to(unit="µs").values[through]
    if not np.isfinite(end).any():
        raise ValueError(f"no neutron reaches {second!r}")
    return {
        "ylim": _with_margin(near, far),
        "xlim": _with_margin(float(np.nanmin(start)), float(np.nanpercentile(end, percentile))),
    }


def _with_margin(low: float, high: float) -> tuple[float, float]:
    """Return ``low`` to ``high`` widened by ``ZOOM_MARGIN`` of their span on either side."""
    margin = ZOOM_MARGIN * (high - low)
    return low - margin, high + margin
