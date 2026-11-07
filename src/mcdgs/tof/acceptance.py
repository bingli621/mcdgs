"""The chopper acceptance diagram of any :class:`.Instrument`.

scippneutron's chopper cascade (:mod:`scippneutron.tof.chopper_cascade`) follows the outline of
a pulse, a polygon in time and wavelength, through the choppers: each one cuts away what reaches
it while it is closed. Moved back to the source, the polygon left after each chopper shows which
times and wavelengths it accepts, and the last one what the whole cascade lets through (J.R.D.
Copley, https://doi.org/10.1016/S0168-9002(03)01731-5).

It takes the same choppers at the same distances as :func:`.build_tof`, through
:func:`.disk_chopper`. The pulse is a rectangle: from 0 to the ESS pulse length, over the
source's band, ``Lmin`` to ``Lmax``, and repeated for the source's ``n_pulses``, unless given.
Unlike tof it has no neutrons, so it is exact and fast, but it does not show the shape of the
pulse.
"""

import math

import matplotlib.pyplot as plt
import numpy as np
import scipp as sc
from matplotlib.figure import Figure
from matplotlib.patches import Polygon
from scippneutron.chopper.disk_chopper import DiskChopper
from scippneutron.tof.chopper_cascade import Chopper, FrameSequence

from mcdgs.core.base import Instrument
from mcdgs.core.components import NXdisk_chopper
from mcdgs.core.components.source import ESS_SOURCE_LENGTH
from mcdgs.core.utils import evaluate, speed_from_wavelength

from .builder import beamline, disk_chopper


def chopper_cascade(
    instrument: Instrument,
    *,
    wavelength_min: float | None = None,
    wavelength_max: float | None = None,
    time_min: float | None = None,
    time_max: float | None = None,
    n_pulses: int | None = None,
) -> tuple[FrameSequence, list[str]]:
    """Follow the pulse of an instrument through its choppers.

    Parameters
    ----------
    instrument:
        The instrument, at its current parameters.
    wavelength_min, wavelength_max:
        Band of the pulse, in Å. ``None`` takes the source's ``Lmin`` or ``Lmax``.
    time_min, time_max:
        Times at the source the pulse starts and ends, in s. ``None`` takes 0 or the ESS pulse
        length.
    n_pulses:
        Number of pulses, one source period apart. ``None`` takes the source's ``n_pulses``, at
        least 1.

    Returns
    -------
    :
        The frames: the pulse at the source, then the one after each chopper, in beam order;
        and the name of each frame, ``"source"`` then the choppers'.

    Raises
    ------
    ValueError
        If the instrument has no ``ESS_butterfly`` or more than one, a chopper is parked, the
        pulse ends before it starts, or ``n_pulses`` is below 1.
    """
    source, after_source = beamline(instrument)
    params = instrument.numeric_instr_params()
    l_min = float(evaluate(source.Lmin, params) if wavelength_min is None else wavelength_min)
    l_max = float(evaluate(source.Lmax, params) if wavelength_max is None else wavelength_max)
    frequency = sc.scalar(float(source.frequency), unit="Hz")
    t_min = 0.0 if time_min is None else float(time_min)
    t_max = ESS_SOURCE_LENGTH if time_max is None else float(time_max)
    if t_max <= t_min:
        raise ValueError(f"the pulse ends at {t_max} s, not after it starts at {t_min} s")
    if n_pulses is None:
        n_pulses = max(source.n_pulses, 1)  # McStas takes 0 as 1
    elif n_pulses < 1:
        raise ValueError(f"n_pulses = {n_pulses}: there must be at least one pulse")

    disks = [
        (name, disk_chopper(component, distance, params), distance)
        for name, component, distance in after_source
        if isinstance(component, NXdisk_chopper)
    ]
    # the slowest neutron of the last pulse passes the last chopper by this time, in s
    last_distance = max((distance for _, _, distance in disks), default=0.0)
    t_end = (n_pulses - 1) / frequency.value + t_max
    t_end += last_distance / speed_from_wavelength(l_max)
    choppers = [_openings(disk, t_min, t_end) for _, disk, _ in disks]

    pulse = FrameSequence.from_source_pulse(
        time_min=sc.scalar(t_min, unit="s"),
        time_max=sc.scalar(t_max, unit="s"),
        wavelength_min=sc.scalar(l_min, unit="angstrom"),
        wavelength_max=sc.scalar(l_max, unit="angstrom"),
        npulses=n_pulses,
        pulse_period=1 / frequency,
    )
    return pulse.chop(choppers), ["source"] + [name for name, _, _ in disks]


def _openings(disk: DiskChopper, t_start: float, t_end: float) -> Chopper:
    """Return the openings of a disk chopper from ``t_start`` to ``t_end``, as a cascade chopper.

    The openings of one rotation, repeated every rotation over the time from ``t_start`` to
    ``t_end``.
    :meth:`DiskChopper.time_offset_open` and :meth:`Chopper.from_disk_chopper` repeat them per
    source period instead, so they need a frequency that is a whole multiple or fraction of the
    source's; the pulse-shaping choppers of T-Rex, at 13.5 times 14 Hz, are not. Their opening
    times also follow the delay, so they need not start within the first period.

    Parameters
    ----------
    disk:
        The scippneutron disk chopper.
    t_start:
        Earliest time a neutron reaches the chopper, in s.
    t_end:
        Latest time a neutron reaches the chopper, in s.

    Returns
    -------
    :
        The chopper, at the distance of the disk's axle.
    """
    period = 1 / abs(disk.frequency.to(unit="Hz").value)  # s, one rotation
    # the edge that passes the beam first opens a slit, as in time_offset_open and _close
    first_edge, last_edge = (
        (disk.slit_begin, disk.slit_end) if disk.is_clockwise else (disk.slit_end, disk.slit_begin)
    )
    # scippneutron gives the rotation before as well, then this one: keep this one
    n_slits = first_edge.sizes["slit"]
    t_open = disk.time_offset_angle_at_beam(angle=first_edge).to(unit="s")["slit", n_slits:]
    t_close = disk.time_offset_angle_at_beam(angle=last_edge).to(unit="s")["slit", n_slits:]
    first = math.floor((t_start - t_close.max().value) / period)
    last = math.ceil((t_end - t_open.min().value) / period)
    shifts = sc.arange("rotation", first, last + 1) * sc.scalar(period, unit="s")
    return Chopper(
        distance=sc.norm(disk.axle_position),
        time_open=(shifts + t_open).flatten(to=t_open.dim),
        time_close=(shifts + t_close).flatten(to=t_close.dim),
    )


def acceptance_diagram(
    instrument: Instrument,
    *,
    wavelength_min: float | None = None,
    wavelength_max: float | None = None,
    time_min: float | None = None,
    time_max: float | None = None,
    n_pulses: int | None = None,
) -> Figure:
    """Draw the chopper acceptance diagram of an instrument.

    The pulse at the source and what each chopper accepts of it, moved back to the source,
    time on x and wavelength on y. Each chopper is drawn over the ones before it, in turbo from the
    source to the last chopper, so the last color is what reaches the end of the cascade. The
    title is the instrument's :meth:`.Instrument.running_conditions`.

    Parameters
    ----------
    instrument:
        The instrument, at its current parameters.
    wavelength_min, wavelength_max:
        Band of the pulse, in Å. ``None`` takes the source's ``Lmin`` or ``Lmax``.
    time_min, time_max:
        Times at the source the pulse starts and ends, in s. ``None`` takes 0 or the ESS pulse
        length.
    n_pulses:
        Number of pulses, one source period apart. ``None`` takes the source's ``n_pulses``, at
        least 1.

    Returns
    -------
    :
        The figure.

    Raises
    ------
    ValueError
        If the instrument has no ``ESS_butterfly`` or more than one, a chopper is parked, the
        pulse ends before it starts, or ``n_pulses`` is below 1.
    """
    frames, names = chopper_cascade(
        instrument,
        wavelength_min=wavelength_min,
        wavelength_max=wavelength_max,
        time_min=time_min,
        time_max=time_max,
        n_pulses=n_pulses,
    )
    at_source = frames.frames[0].distance
    colors = plt.get_cmap("turbo")(np.linspace(0, 1, len(frames.frames)))

    fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
    for frame, name, color in zip(frames.frames, names, colors, strict=True):
        label = name if name == "source" else f"{name} ({frame.distance.value:.2f} m)"
        ax.plot([], [], color=color, linewidth=6, label=label)  # one legend entry per frame
        for subframe in frame.propagate_to(at_source).subframes:
            wavelength = subframe.wavelength.to(unit="angstrom").values
            time = subframe.time.to(unit="ms").values
            ax.add_patch(Polygon(np.stack((time, wavelength), axis=1), closed=True, color=color))
    ax.autoscale_view()
    ax.set_xlabel("Time at the source (ms)")
    ax.set_ylabel("Wavelength (Å)")
    ax.set_title(instrument.running_conditions())
    ax.grid(True, linewidth=0.5, alpha=0.5)
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)
    return fig
