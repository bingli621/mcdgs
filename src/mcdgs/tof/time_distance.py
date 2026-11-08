"""The time-distance diagram of any :class:`.Instrument`, from its tof model.

:func:`time_distance_diagram` builds the tof model of the instrument (:func:`.build_tof`), runs
it and draws tof's time-distance diagram: each neutron is a ray from the moderator, through the
openings of the choppers, to the detectors, colored by its wavelength, and the neutrons a
chopper stops are drawn in gray up to it. Unlike the acceptance diagram (:mod:`.acceptance`), it
follows neutrons, so it shows the shape of the pulse, with Monte Carlo noise.
"""

from matplotlib.figure import Figure

from mcdgs.core.base import Instrument
from mcdgs.tof.builder import build_tof


def time_distance_diagram(
    instrument: Instrument,
    neutrons: int = 1_000_000,
    *,
    seed: int | None = None,
    n_pulses: int | None = None,
    visible_rays: int = 2000,
    blocked_rays: int = 5000,
    figsize: tuple[float, float] = (10, 6),
    cmap: str = "turbo",
) -> Figure:
    """Build the tof model of an instrument, run it and draw its time-distance diagram.

    The title is the instrument's :meth:`.Instrument.running_conditions`. tof draws the time in
    µs.

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
        Number of pulses, one source period apart. ``None`` takes the source's ``n_pulses``, at
        least 1.
    visible_rays:
        Most rays drawn, per pulse, of the neutrons that reach the last detector.
    blocked_rays:
        Most rays drawn in gray, per pulse, of the neutrons a chopper stops. 0 draws none.
    figsize:
        Size of the figure, in inches.
    cmap:
        Colormap of the wavelength.

    Returns
    -------
    :
        The figure.

    Raises
    ------
    ValueError
        If the instrument has no ``ESS_butterfly`` or more than one, or a chopper is parked.
    """
    result = build_tof(instrument, neutrons=neutrons, seed=seed, n_pulses=n_pulses).run()
    plot = result.plot(
        visible_rays=visible_rays,
        blocked_rays=blocked_rays,
        figsize=figsize,
        cmap=cmap,
        seed=seed,
        title=instrument.running_conditions(),
    )
    plot.fig.tight_layout()
    return plot.fig
