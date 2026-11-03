"""Build the tof model of any :class:`.Instrument`.

tof follows each neutron along a single flight path from the source, through the choppers, to
the detectors. It has no optics: guides, benders and slits transmit everything. So the builder
keeps only the components tof knows, by type:

- the ``ESS_butterfly`` becomes the :class:`tof.Source`, which must be the only source;
- each ``NXdisk_chopper`` becomes a :class:`tof.Chopper`;
- each ``Monitor_nD``, and each component named in ``detectors`` (such as the ``Arm`` at the
  sample position), becomes a :class:`tof.Detector`, which records the neutrons reaching it.

The optics and the markers (``SKIPPED``) are left out. Any other component raises, so that a new
component class is never ignored without a decision.

Distances: a component's distance is the straight line from the source to its position
(:meth:`.Instrument.positions`). For T-Rex this is within 1 mm of the path along the curved
guide, under 1 % of the shortest chopper opening at 2 Å. The tof source is put at 0 m, where
McStas has the moderator, rather than at tof's default for the ESS pulse.
"""

from collections.abc import Callable, Collection, Mapping

import numpy as np
import scipp as sc
import tof

from mcdgs.core.base import Instrument
from mcdgs.core.components import (
    Arm,
    Component,
    ESS_butterfly,
    Guide_gravity,
    Monitor_nD,
    NXdisk_chopper,
    Pol_bender,
    Progress_bar,
    Slit,
)
from mcdgs.core.utils import evaluate

# what tof leaves out: the optics, which it takes as transmitting everything, and the markers
SKIPPED: frozenset[type[Component]] = frozenset(
    {Arm, Guide_gravity, Pol_bender, Progress_bar, Slit}
)

# a converter takes the component, its distance from the source in m, and the instrument
# parameter values; it returns the tof component
Converter = Callable[[Component, float, Mapping[str, float]], tof.Component]


def _numeric_params(instrument: Instrument) -> dict[str, float]:
    """Return the instrument parameters that are numbers, by name.

    Only numbers can stand in a component's values; a string parameter cannot.

    Parameters
    ----------
    instrument:
        The instrument.

    Returns
    -------
    :
        Value by McStas name.
    """
    return {p.name: p.value for p in instrument.instr_params if isinstance(p.value, int | float)}


def _source(
    component: ESS_butterfly,
    params: Mapping[str, float],
    neutrons: int,
    seed: int | None,
) -> tof.Source:
    """Return the tof source of the ESS butterfly moderator.

    tof samples its standard ESS pulse, the same for every beamline: ``sector``, ``beamline``,
    ``cold_frac``, ``yheight``, the focusing and the performance factors are not used.

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

    Returns
    -------
    :
        The source at 0 m, from ``Lmin`` to ``Lmax``, with ``n_pulses`` pulses (at least 1) at
        the moderator's frequency.
    """
    return tof.Source(
        facility="ess",
        neutrons=neutrons,
        pulses=max(component.n_pulses, 1),  # McStas takes 0 as 1
        frequency=sc.scalar(float(component.frequency), unit="Hz"),
        wmin=sc.scalar(float(evaluate(component.Lmin, params)), unit="angstrom"),
        wmax=sc.scalar(float(evaluate(component.Lmax, params)), unit="angstrom"),
        seed=seed,
        distance=sc.scalar(0.0, unit="m"),
    )


def _chopper(
    component: NXdisk_chopper, distance: float, params: Mapping[str, float]
) -> tof.Chopper:
    """Return the tof chopper of a disk chopper.

    Parameters
    ----------
    component:
        The disk chopper.
    distance:
        Its distance from the source, in m.
    params:
        Instrument parameter values, for ``nu``, ``delay`` and the angles given as expressions.

    Returns
    -------
    :
        The tof chopper.

    Raises
    ------
    NotImplementedError
        Always, until the conversion of the phase from McStas to tof is written.
    """
    raise NotImplementedError(f"{component.name}: NXdisk_chopper is not yet converted to tof")


def _detector(component: Component, distance: float, params: Mapping[str, float]) -> tof.Detector:
    """Return a tof detector at the place of a component.

    Parameters
    ----------
    component:
        The component recording the neutrons, such as a monitor or the sample position.
    distance:
        Its distance from the source, in m.
    params:
        Instrument parameter values (not used: a detector is only a distance).

    Returns
    -------
    :
        The detector, named after the component.
    """
    return tof.Detector(distance=sc.scalar(distance, unit="m"), name=component.name)


CONVERTERS: dict[type[Component], Converter] = {
    NXdisk_chopper: _chopper,
    Monitor_nD: _detector,
}


def build_tof(
    instrument: Instrument,
    *,
    neutrons: int = 1_000_000,
    detectors: Collection[str] = (),
    seed: int | None = None,
) -> tof.Model:
    """Build the tof model of an instrument.

    Parameters
    ----------
    instrument:
        The instrument to build, at its current parameters.
    neutrons:
        Number of neutrons per pulse.
    detectors:
        Names of components to record the neutrons at, in addition to every ``Monitor_nD``,
        such as ``("sample_pos",)``. Each must be after the source.
    seed:
        Seed of the random numbers of the source. If ``None``, every build differs.

    Returns
    -------
    :
        The tof model: the source, then the choppers and detectors in beam order.

    Raises
    ------
    ValueError
        If the instrument has no ``ESS_butterfly`` or more than one, if a name in ``detectors``
        is not a component after the source, or if a component is neither converted nor
        skipped.
    """
    components = instrument.components
    sources = [name for name, c in components.items() if isinstance(c, ESS_butterfly)]
    if len(sources) != 1:
        raise ValueError(f"{instrument.name}: tof needs one ESS_butterfly, found {sources}")
    source_name = sources[0]
    after_source = list(components)[list(components).index(source_name) + 1 :]
    if unknown := sorted(set(detectors) - set(after_source)):
        raise ValueError(f"{instrument.name}: no component {unknown} after the source")

    params = _numeric_params(instrument)
    positions = {name: np.array(xyz) for name, xyz in instrument.positions().items()}
    model = tof.Model(
        source=_source(components[source_name], params, neutrons, seed),
    )
    for name in after_source:
        component = components[name]
        converter = _detector if name in detectors else CONVERTERS.get(type(component))
        if converter is None:
            if type(component) in SKIPPED:
                continue
            raise ValueError(
                f"{name} ({component.component_name}): no tof conversion; add one to "
                "CONVERTERS or the class to SKIPPED"
            )
        distance = float(np.linalg.norm(positions[name] - positions[source_name]))
        model.add(converter(component, distance, params))
    return model
