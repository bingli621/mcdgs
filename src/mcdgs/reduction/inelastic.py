"""The energy and momentum transfer of neutrons scattered by the sample, in direct geometry.

``QE_GRAPH`` is a scipp graph built from the kernels of scippneutron. It takes the incident
energy of a sub-pulse and when it reaches the sample, and where and when a neutron reaches the
detector, and gives the final energy, the energy transfer and the momentum transfer of the
neutron. The energy transfer is ``Ei - Ef``, positive when the neutron loses energy, and the
momentum transfer ``ki - kf``.

A pulse of a repetition-rate multiplication instrument holds several sub-pulses, each with its
own incident energy. :func:`arrival_window` gives, per pixel and sub-pulse, when the neutrons of
an energy-transfer range reach the detector, and :func:`assign_subpulse` gives each neutron the
sub-pulse of the window it arrives in.
"""

from typing import cast

import numpy as np
import scipp as sc
from scipp.coords.graph import GraphDict
from scipp.coords.rule import Kernel
from scippneutron.conversion import beamline as beamline_kernels
from scippneutron.conversion import tof as tof_kernels

__all__ = [
    "QE_GRAPH",
    "arrival_window",
    "assign_subpulse",
    "final_wavelength",
    "time_at_sample",
]

DEFAULT_ENERGY_TRANSFER_RATIO = (-0.8, 0.8)  # E / Ei, from energy gain to energy loss


def time_at_sample(
    time_at_monitor: sc.Variable,
    monitor_position: sc.Variable,
    sample_position: sc.Variable,
    incident_energy: sc.Variable,
) -> sc.Variable:
    """Return when a neutron reaches the sample, from when it passes a monitor before it.

    Parameters
    ----------
    time_at_monitor:
        When the neutron passes the monitor.
    monitor_position, sample_position:
        Where the monitor and the sample are, as vectors.
    incident_energy:
        The energy of the neutron.

    Returns
    -------
    :
        The time at the sample, in the unit of ``time_at_monitor``.

    Examples
    --------
    A neutron of 2 Å flies at 1978 m/s, so 1 m takes 0.50556 ms:

    >>> t = time_at_sample(
    ...     sc.scalar(100.0, unit="ms"),
    ...     sc.vector([0.0, 0.0, 0.0], unit="m"),
    ...     sc.vector([0.0, 0.0, 1.0], unit="m"),
    ...     tof_kernels.energy_from_wavelength(wavelength=sc.scalar(2.0, unit="Å")),
    ... )
    >>> print(f"{t.value:.5f} {t.unit}")
    100.50556 ms
    """
    wavelength = tof_kernels.wavelength_from_energy(energy=incident_energy)
    distance = sc.norm(sample_position - monitor_position)
    flight = tof_kernels.tof_from_wavelength(wavelength=wavelength, Ltotal=distance)
    return time_at_monitor + flight.to(unit=time_at_monitor.unit)


# --8<-- [start:qe_graph]
def final_wavelength(
    time_at_sample: sc.Variable, time_at_detector: sc.Variable, L2: sc.Variable
) -> sc.Variable:
    """Return the wavelength of a neutron from its flight from the sample to the detector."""
    return tof_kernels.wavelength_from_tof(tof=time_at_detector - time_at_sample, Ltotal=L2)


def _incident_wavelength(incident_energy: sc.Variable) -> sc.Variable:
    """Return the wavelength of the incident neutron."""
    return tof_kernels.wavelength_from_energy(energy=incident_energy)


def _incident_wavevector(
    incident_wavelength: sc.Variable, incident_beam: sc.Variable
) -> sc.Variable:
    """Return the wavevector ki of the incident neutron, along the incident beam."""
    return tof_kernels.wavevector_from_wavelength(
        wavelength=incident_wavelength, beam=incident_beam
    )


def _final_wavevector(final_wavelength: sc.Variable, scattered_beam: sc.Variable) -> sc.Variable:
    """Return the wavevector kf of the scattered neutron, along the scattered beam."""
    return tof_kernels.wavevector_from_wavelength(wavelength=final_wavelength, beam=scattered_beam)


def _momentum_transfer(
    incident_wavevector: sc.Variable, final_wavevector: sc.Variable
) -> sc.Variable:
    """Return the momentum transfer ``ki - kf``."""
    return incident_wavevector - final_wavevector.to(unit=incident_wavevector.unit)


def _momentum_transfer_magnitude(Q_vec: sc.Variable) -> sc.Variable:
    """Return the magnitude of the momentum transfer."""
    return sc.norm(Q_vec)


def _momentum_transfer_components(Q_vec: sc.Variable) -> dict[str, sc.Variable]:
    """Return the x, y and z components of the momentum transfer, in the lab frame."""
    return {"Qx": Q_vec.fields.x.copy(), "Qy": Q_vec.fields.y.copy(), "Qz": Q_vec.fields.z.copy()}


# scippneutron types them for data arrays and datasets too, scipp's graph only for variables
_incident_beam = cast(Kernel, beamline_kernels.straight_incident_beam)
_scattered_beam = cast(Kernel, beamline_kernels.straight_scattered_beam)
_L2 = cast(Kernel, beamline_kernels.L2)

# from the incident energy of a sub-pulse and its time at the sample, and the position and time
# of a neutron at the detector, to its energy and momentum transfer
QE_GRAPH: GraphDict = {
    "incident_beam": _incident_beam,  # from source_position to sample_position
    "scattered_beam": _scattered_beam,  # from sample_position to position
    "L2": _L2,
    "final_wavelength": final_wavelength,
    "final_energy": tof_kernels.final_energy_from_wavelength,
    "energy_transfer": tof_kernels.energy_transfer_from_energies,  # Ei - Ef
    "incident_wavelength": _incident_wavelength,
    "incident_wavevector": _incident_wavevector,
    "final_wavevector": _final_wavevector,
    "Q_vec": _momentum_transfer,  # ki - kf
    "Q": _momentum_transfer_magnitude,
    ("Qx", "Qy", "Qz"): _momentum_transfer_components,
}
# --8<-- [end:qe_graph]


def arrival_window(
    time_at_sample: sc.Variable,
    L2: sc.Variable,
    incident_energy: sc.Variable,
    energy_transfer_ratio: tuple[float, float] = DEFAULT_ENERGY_TRANSFER_RATIO,
) -> dict[str, sc.Variable]:
    """Return when the neutrons of a range of energy transfer reach the detector.

    The range is that of ``E / Ei``, with ``E = Ei - Ef``: a neutron that gains energy, with a
    negative ratio, arrives first, and one that loses energy last.

    Parameters
    ----------
    time_at_sample:
        When each sub-pulse reaches the sample.
    L2:
        The distance from the sample to each pixel.
    incident_energy:
        The incident energy of each sub-pulse, with the dimensions of ``time_at_sample``.
    energy_transfer_ratio:
        The smallest and the largest ``E / Ei``. The largest must be below 1: a neutron that
        loses all its energy never arrives.

    Returns
    -------
    :
        ``earliest_arrival`` and ``latest_arrival``, by name, with the dimensions of
        ``time_at_sample`` and ``L2``, in the unit of ``time_at_sample``.

    Raises
    ------
    ValueError
        If the range is empty or its end is not below 1.

    Examples
    --------
    A neutron of 2 Å flies 3 m in 1.5167 ms; with E / Ei = 0.75 it has half the speed, and takes
    twice as long:

    >>> window = arrival_window(
    ...     sc.scalar(0.0, unit="ms"),
    ...     sc.scalar(3.0, unit="m"),
    ...     tof_kernels.energy_from_wavelength(wavelength=sc.scalar(2.0, unit="Å")),
    ...     energy_transfer_ratio=(0.0, 0.75),
    ... )
    >>> print(f"{window['earliest_arrival'].value:.4f}, {window['latest_arrival'].value:.4f}")
    1.5167, 3.0333
    """
    smallest, largest = energy_transfer_ratio
    if not smallest < largest < 1:
        raise ValueError(
            f"energy_transfer_ratio must rise and end below 1, got {energy_transfer_ratio}"
        )
    window = {}
    for name, ratio in (("earliest_arrival", smallest), ("latest_arrival", largest)):
        final_energy = incident_energy * (1 - ratio)
        wavelength = tof_kernels.wavelength_from_energy(energy=final_energy)
        flight = tof_kernels.tof_from_wavelength(wavelength=wavelength, Ltotal=L2)
        window[name] = time_at_sample + flight.to(unit=time_at_sample.unit)
    return window


def assign_subpulse(
    pixel_id: sc.Variable,
    time_at_detector: sc.Variable,
    window_pixel_id: sc.Variable,
    earliest_arrival: sc.Variable,
    latest_arrival: sc.Variable,
) -> dict[str, sc.Variable]:
    """Give each neutron the sub-pulse of the arrival window it falls in.

    A neutron belongs to the last sub-pulse whose window has opened when it arrives, at its
    pixel; it is inside if that window has not closed yet. Where the windows of two sub-pulses
    overlap, the later one wins. A neutron that arrives before the first window opens belongs to
    the first sub-pulse and is outside.

    The times are those within one pulse: frames are not unwrapped.

    Parameters
    ----------
    pixel_id, time_at_detector:
        The pixel and the time of each neutron, along one dimension.
    window_pixel_id:
        The pixels of the windows, along their pixel dimension. Every pixel of ``pixel_id`` must
        be in it.
    earliest_arrival, latest_arrival:
        The windows of :func:`arrival_window`, along the pixel dimension of ``window_pixel_id``
        and the sub-pulse dimension, the sub-pulses in time order.

    Returns
    -------
    :
        ``subpulse``, the index of the sub-pulse of each neutron, and ``inside``, whether it falls
        in the window of that sub-pulse, by name, along the dimension of ``pixel_id``.

    Raises
    ------
    ValueError
        If the windows do not have two dimensions, or a pixel of a neutron has no window.
    """
    (pixel_dim,) = window_pixel_id.dims
    if earliest_arrival.ndim != 2 or pixel_dim not in earliest_arrival.dims:
        raise ValueError(
            f"The windows must have the dimension {pixel_dim!r} and a sub-pulse one, "
            f"got {earliest_arrival.dims}"
        )
    unit = earliest_arrival.unit
    dims = [pixel_dim, *(d for d in earliest_arrival.dims if d != pixel_dim)]
    # the windows of each neutron's pixel: (neutron, sub-pulse)
    order = np.argsort(window_pixel_id.values)
    sorted_ids = window_pixel_id.values[order]
    found = np.searchsorted(sorted_ids, pixel_id.values).clip(max=sorted_ids.size - 1)
    if not np.array_equal(sorted_ids[found], pixel_id.values):
        raise ValueError("A pixel of a neutron has no arrival window")
    rows = order[found]
    opens = earliest_arrival.transpose(dims).values[rows]
    closes = latest_arrival.to(unit=unit).transpose(dims).values[rows]
    times = time_at_detector.to(unit=unit).values[:, None]

    subpulse = ((times >= opens).sum(axis=1) - 1).clip(min=0)
    neutron = np.arange(subpulse.size)
    inside = (times[:, 0] >= opens[neutron, subpulse]) & (times[:, 0] < closes[neutron, subpulse])
    (dim,) = pixel_id.dims
    return {
        "subpulse": sc.array(dims=[dim], values=subpulse),
        "inside": sc.array(dims=[dim], values=inside),
    }
