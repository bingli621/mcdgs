"""Build the tof model of any :class:`.Instrument`, and give its disks as scippneutron choppers.

tof follows each neutron along a single flight path from the source, through the choppers, to
the detectors. It has no optics: guides, benders and slits transmit everything. So the builder
takes only the components tof understands, by type:

- the ``ESS_butterfly``, or the ``Polygon_ESS_butterfly`` that extends it, becomes the
  :class:`tof.Source`, which must be the only source;
- each ``NXdisk_chopper`` becomes a :class:`tof.Chopper`, through scippneutron's
  :class:`DiskChopper` (:func:`disk_chopper`), as both follow NeXus;
- each ``Monitor_nD`` becomes a :class:`tof.Detector`, which records the neutrons reaching it.
  An instrument without one gives a model without detectors. A curved one, such as a banana,
  is centered on the sample: its detector is its radius further on.

Every other component, such as the optics, the markers or a sample, is left out.

:func:`get_choppers` gives the scippneutron :class:`DiskChopper` of each disk chopper by name,
which scippneutron can draw and :mod:`.acceptance` follows the pulse through.
:func:`beamline` gives the source and the distance of each component after it, for both.

The source: tof has a pulse of the ESS moderator for each instrument, sampled at its beam port.
The builder takes the one of the moderator's beam port, named by its ``sector`` and
``beamline`` (``"W7"``); a beam port it does not know gets tof's standard ESS pulse, ``"ess"``.

Distances: a component's distance is the straight line from the source to its position
(:meth:`.Instrument.positions`). For T-Rex this is within 1 mm of the path along the curved
guide, under 1 % of the shortest chopper opening at 2 Å. McStas has the moderator at 0 m, and
the tof source stays where tof sampled its pulse, 5 cm from it: its times are those of the
neutrons there, so putting it at 0 m would add the flight over those 5 cm twice.
"""

from collections.abc import Mapping

import numpy as np
import scipp as sc
import tof
from scippneutron.chopper.disk_chopper import DiskChopper

from mcdgs.core.base import Instrument
from mcdgs.core.components import Component, ESS_butterfly, Monitor_nD, NXdisk_chopper
from mcdgs.core.utils import evaluate

# tof's pulse of the ESS moderator by beam port, ``sector`` and ``beamline``; else ``"ess"``
TOF_SOURCES = {"W7": "ess-trex", "W3": "ess-cspec"}


def source_pulses(source: ESS_butterfly, n_pulses: int | None = None) -> int:
    """Return the number of pulses of a source, as McStas takes it.

    Parameters
    ----------
    source:
        The moderator.
    n_pulses:
        The number of pulses. If ``None``, the moderator's ``n_pulses``, at least 1: McStas
        takes 0 as 1.

    Returns
    -------
    :
        The number of pulses.

    Raises
    ------
    ValueError
        If ``n_pulses`` is below 1.
    """
    if n_pulses is None:
        return max(source.n_pulses, 1)
    if n_pulses < 1:
        raise ValueError(f"n_pulses = {n_pulses}: there must be at least one pulse")
    return n_pulses


def disk_chopper(
    component: NXdisk_chopper, distance: float, params: Mapping[str, float]
) -> DiskChopper:
    """Return the scippneutron disk chopper of a McStas disk chopper.

    Both follow NeXus ``NXdisk_chopper``. On the beam axis, McStas has the angle under the beam,
    from the zero mark in the sense of ``slit_edges``, at
    ``beam_angle - 360 * nu * (t - delay)``; scippneutron at
    ``beam_position + phase - 360 * frequency * t``. So the slit edges, the signed frequency and
    the beam angle carry over, and the phase is ``360 * nu * delay``, as
    :meth:`DiskChopper.from_nexus` computes it from a delay.

    Only the timing on the beam axis counts for tof. ``zero_angle`` only says where the spindle
    lies, so it does not change it; ``xwidth``, ``jitter`` and ``abs_out`` are left out. The
    ``radius`` and, as the height of every slit, ``yheight`` carry over so that scippneutron
    draws the disk to scale, but tof does not use them; a ``yheight`` of 0 reaches the spindle,
    so the slits are as high as the radius.

    Parameters
    ----------
    component:
        The McStas disk chopper.
    distance:
        Its distance from the source, in m: the ``z`` of the axle.
    params:
        Instrument parameter values, for ``nu``, ``delay`` and ``beam_angle`` given as
        expressions.

    Returns
    -------
    :
        The disk chopper.

    Raises
    ------
    ValueError
        If ``nu`` is 0. McStas parks the disk, which tof cannot: it takes a disk at rest as open.

    Examples
    --------
    >>> chopper = NXdisk_chopper(name="c", slit_edges=[-10, 10], nu="f", delay=0.01)
    >>> disk = disk_chopper(chopper, 20.0, {"f": 14})
    >>> float(disk.frequency.value), float(disk.phase.value), disk.slit_end.values.tolist()
    (14.0, 50.4, [10.0])
    """
    nu = float(evaluate(component.nu, params))
    if nu == 0:
        raise ValueError(f"{component.name}: a parked disk (nu = 0) cannot be modeled in tof")
    delay = float(evaluate(component.delay, params))
    edges = np.array(component.slit_edges, float)
    n_slits = len(edges) // 2
    slit_height = component.yheight if component.yheight > 0 else component.radius
    return DiskChopper(
        axle_position=sc.vector([0.0, 0.0, distance], unit="m"),
        frequency=sc.scalar(nu, unit="Hz"),
        beam_position=sc.scalar(float(evaluate(component.beam_angle, params)), unit="deg"),
        phase=sc.scalar(360 * nu * delay, unit="deg"),
        slit_begin=sc.array(dims=["slit"], values=edges[0::2], unit="deg"),
        slit_end=sc.array(dims=["slit"], values=edges[1::2], unit="deg"),
        slit_height=sc.full(dims=["slit"], shape=[n_slits], value=float(slit_height), unit="m"),
        radius=sc.scalar(float(component.radius), unit="m"),
    )


def get_choppers(instrument: Instrument) -> dict[str, DiskChopper]:
    """Return the scippneutron disk choppers of an instrument, through :func:`disk_chopper`.

    Parameters
    ----------
    instrument:
        The instrument, at its current parameters.

    Returns
    -------
    :
        The disk chopper of each ``NXdisk_chopper`` after the source, by name, in beam order.

    Raises
    ------
    ValueError
        If the instrument has no ``ESS_butterfly`` or ``Polygon_ESS_butterfly``, or more than
        one; or if a chopper is parked.
    """
    _, after_source = beamline(instrument)
    params = instrument.numeric_instr_params()
    return {
        name: disk_chopper(component, distance, params)
        for name, component, distance in after_source
        if isinstance(component, NXdisk_chopper)
    }


def _source(
    component: ESS_butterfly,
    params: Mapping[str, float],
    neutrons: int,
    seed: int | None,
    n_pulses: int | None = None,
) -> tof.Source:
    """Return the tof source of the ESS butterfly moderator.

    tof samples the pulse of the moderator's beam port, ``sector`` and ``beamline`` such as
    ``"W7"``, if :data:`TOF_SOURCES` has it, else its standard ESS pulse, at the distance tof
    sampled it at. ``cold_frac``, ``yheight``, the focusing and the performance factors are
    not used.

    Parameters
    ----------
    component:
        The moderator.
    params:
        Instrument parameter values, for ``Lmin`` and ``Lmax`` given as expressions.
    neutrons:
        Number of neutrons per pulse.
    seed:
        Seed of the random numbers. If ``None``, every run differs.
    n_pulses:
        Number of pulses, as :func:`source_pulses`.

    Returns
    -------
    :
        The source, from ``Lmin`` to ``Lmax``, with ``n_pulses`` pulses at the moderator's
        frequency.
    """
    return tof.Source(
        facility=TOF_SOURCES.get(f"{component.sector}{component.beamline}", "ess"),
        neutrons=neutrons,
        pulses=source_pulses(component, n_pulses),
        frequency=sc.scalar(float(component.frequency), unit="Hz"),
        wmin=sc.scalar(float(evaluate(component.Lmin, params)), unit="angstrom"),
        wmax=sc.scalar(float(evaluate(component.Lmax, params)), unit="angstrom"),
        seed=seed,
    )


def build_tof(
    instrument: Instrument,
    *,
    neutrons: int = 1_000_000,
    seed: int | None = None,
    n_pulses: int | None = None,
) -> tof.Model:
    """Build the tof model of an instrument.

    Parameters
    ----------
    instrument:
        The instrument to build, at its current parameters.
    neutrons:
        Number of neutrons per pulse.
    seed:
        Seed of the random numbers of the source. If ``None``, every build differs.
    n_pulses:
        Number of pulses of the source, as :func:`source_pulses`.

    Returns
    -------
    :
        The tof model: the source, then the choppers and the detectors (one per ``Monitor_nD``)
        in beam order, each named after its component.

    Raises
    ------
    ValueError
        If the instrument has no ``ESS_butterfly`` or ``Polygon_ESS_butterfly``, or more than
        one; if a chopper is parked; or if ``n_pulses`` is below 1.
    """
    source, after_source = beamline(instrument)
    params = instrument.numeric_instr_params()
    model = tof.Model(source=_source(source, params, neutrons, seed, n_pulses))
    for name, component, distance in after_source:
        match component:
            case NXdisk_chopper():
                # tof's conversion turns the signed frequency into a direction, and the phase
                disk = disk_chopper(component, distance, params)
                model.add(tof.Chopper.from_diskchopper(disk, name=name))
            case Monitor_nD():  # a curved one is centered on the sample: its radius further on
                reached = distance + (component.curved_radius or 0)
                model.add(tof.Detector(distance=sc.scalar(reached, unit="m"), name=name))
            case _:
                pass  # not something tof understands
    return model


def beamline(
    instrument: Instrument,
) -> tuple[ESS_butterfly, list[tuple[str, Component, float]]]:
    """Return the source of an instrument and the components after it, with their distances.

    Parameters
    ----------
    instrument:
        The instrument.

    Returns
    -------
    :
        The source, and each component after it in beam order: its name, itself and
        its distance from the source, in m.

    Raises
    ------
    ValueError
        If the instrument has no ``ESS_butterfly`` or ``Polygon_ESS_butterfly``, or more than
        one.
    """
    components = instrument.components
    sources = {name: c for name, c in components.items() if isinstance(c, ESS_butterfly)}
    if len(sources) != 1:
        raise ValueError(
            f"{instrument.name}: tof needs one ESS_butterfly or Polygon_ESS_butterfly, "
            f"found {list(sources)}"
        )
    ((source_name, source),) = sources.items()
    after_source = list(components)[list(components).index(source_name) + 1 :]

    positions = {name: np.array(xyz) for name, xyz in instrument.positions().items()}
    origin = positions[source_name]
    return source, [
        (name, components[name], float(np.linalg.norm(positions[name] - origin)))
        for name in after_source
    ]
