"""The analytic elastic resolution of T-Rex: the width of the elastic line, without McStas.

The width of the elastic line of vanadium comes from three places: how long the pulse-shaping
choppers are open, how long the monochromatic choppers are open, and where in the sample a neutron
scatters. This module computes each from the geometry of :class:`.TRex`, convolves them, and fits
a Gaussian to the result as :func:`.elastic_resolution_incoherent.fit_elastic_lines` fits the
simulated lines. It needs no McStas and no random numbers.

:func:`literature_resolution` gives, for comparison, the same width from eq. (3) of Lechner
(2000), the main reference here, as written on the elastic line: the opening time of each pair
from its slit alone, ``a / (2 u)`` as in eq. (1) of Gaspar (2007), and the three terms added in
quadrature. The model keeps eq. (3), and differs from it in two places, step 2 and step 5 below.
It goes in six steps.

**1. The flight of an elastic neutron.** A neutron passes the pulse-shaping pair at ``t_P`` and
the monochromatic pair at ``t_M``, ``L_PM`` apart, which fixes its speed,
``v = L_PM / (t_M - t_P)``. It reaches the sample, ``L_MS`` after the monochromatic pair, at
``t_S = t_M + L_MS / v``, and the detector, ``L2`` after the sample, at ``t_D = t_S + L2 / v``.
The reduction does not know ``t_P`` and ``t_M`` of each neutron: it gives every neutron of a
sub-pulse the same incident energy ``Ei`` and the same time at the sample ``t_S0``, and takes
``Ef = m (L2 / (t_D - t_S0))² / 2`` and ``E = Ei - Ef``. To first order in the offsets
``δt_P``, ``δt_M`` of the neutron from the centers of the two openings, the flight time it
measures is off by

    δ(t_D - t_S0) = A δt_M - B δt_P,    A = 1 + L_MD / L_PM,    B = L_MD / L_PM,

with ``L_MD = L_MS + L2``, the path from the monochromatic choppers to the detector. Since
``E = Ei - Ef`` and ``dEf / dt = -2 Ef / t``, with ``t = L2 / v0`` and ``Ef = Ei`` on the elastic
line, the energy transfer it computes is off by

    δE = (2 Ei v0 / L2) (A δt_M - B δt_P),

with ``v0`` the speed of the sub-pulse. This is eq. (3) of Lechner (2000) at ``λ = λ0``, with his
``L12`` for ``L_PM`` and his ``L2S + LSD`` for ``L_MD``, and eq. (19) of Gaspar (2007) at
``Ei = Ef``: there
``L0 = L_PM``, and the path from the monochromatic choppers to the detector is ``L2 + L1`` in its
notation. The pulse-shaping choppers weigh in only by ``B``, 0.09 for T-Rex, but they are open
about ten times as long as the monochromatic ones, so their term is as large. Both terms are in
proportion to ``Ei v0``, so FWHM / Ei grows as the square root of Ei.

**2. How long a counter-rotating pair is open.** Each pair is two disks with the same slit,
moving across the beam in opposite directions, each at its rim speed ``u``; :func:`read_pair`
checks this. A neutron at ``x`` across the beam passes if it is inside both slits. The two slits
of arc ``a`` overlap on a window ``a - 2 u |t|`` wide, which shrinks to nothing in ``a / (2 u)``.
If the beam is at least as wide as the slit, the whole window is on the beam, and the pair is open
as a triangle in time with a FWHM of ``a / (2 u)``: eq. (1) of Gaspar (2007),
``Δt = b / (2π f D)`` with ``u = π f D``. A beam narrower than the slit cuts the top off the
triangle: the pair is open as a trapezoid, and longer, up to ``a / u``, a single disk, for a beam
of no width. T-Rex has both: at the pulse-shaping choppers the beam, 60 mm, is narrower than the
slits, 122 and 214 mm, so their openings are long trapezoids; at the monochromatic choppers it,
22.4 mm, is wider than the narrow slit, 15.3 mm, a triangle that also cuts the beam, and narrower
than the wide one, 26.9 mm. :func:`pair_transmission` computes the opening exactly, as the length
of the beam inside both slits at each time, for a beam that fills its window evenly. It also takes
in that the second disk is ``spacing`` further down the beam: its delay is set for the speed of
``lambda0``, so for a sub-pulse of another speed the two slits are not centered at the same
time, by ``spacing (1 / v - 1 / v_lambda0)``. Lechner (2000) gives the opening times without a
formula, but notes that the window of the chopper and the aperture of the guide set them together
(his Sec. 5.2, and Fig. 2); Gaspar (2007) takes eq. (1), the triangle, for every pair.

**3. Where the neutron scatters.** A neutron that scatters at ``(x, y)`` in the sample, ``x``
along the beam and ``y`` across it in the horizontal plane, from its center, into the scattering
angle 2θ flies further by

    δpath = x (1 - cos 2θ) - y sin 2θ,

which the reduction, which puts every neutron at the center, reads as a time ``δpath / v0`` and an
error ``δE = 2 Ei δpath / L2``. This is the ``ΔL`` term of eq. (3) of Lechner (2000), as eq. (20)
of Gaspar (2007) gives it for a rectangle, at ``Ei = Ef``. For the cylinder of radius ``r`` of T-Rex, :func:`sample_path_density` uses that
``δpath`` is the projection of a uniform disk onto a vector of length ``2 sin θ``: a semicircle of
radius ``2 r |sin θ|``. It averages it over the 2θ of the detector, evenly per degree: vanadium
scatters the same into every direction, and the pixels of the banana are of equal angle. Out of the
horizontal plane, the detector spans 25° below to 15° above the beam; that is left out.

**4. The detector.** The banana of the simulation has no depth and lists the exact time of each
neutron, so it adds nothing. A real detector adds its thickness to ``ΔL`` in eq. (3) of Lechner
(2000), as his Table 1 does for 15 mm ³He tubes, or as eq. (9) of Gaspar (2007).

**5. Together.** The three offsets are independent: the moderator pulse, 2.86 ms long, is flat
across the openings, so ``δt_P`` and ``δt_M`` are not correlated, nor is either with where the
neutron scatters. :func:`elastic_resolution` therefore convolves the three distributions of
``δE``. It does not add their widths in quadrature, as eq. (3) of Lechner (2000) does, and
Gaspar (2007) for Gaussians (his footnote 1): the openings are triangles and trapezoids, for which
quadrature gives a line about 25 % too wide. Lechner notes that for pulses that are not Gaussian,
standard deviations may serve better than the FWHM (his Sec. 3).

**6. The width.** As the reduction, it fits a Gaussian to the bins of at least ``fit_fraction``
of the highest, and gives ``FWHM_PER_SIGMA`` times its σ.

What it leaves out: terms of second order in the offsets, the spread of path lengths in the guide
between the choppers, gravity, and the bandwidth choppers. The last cut the first and the last
sub-pulse of a frame, which makes them narrower than the model.

References
----------
R. E. Lechner, "Quasielastic high-resolution time-of-flight spectrometers employing multi-disk
chopper cascades for spallation sources", Proceedings of ICANS-XV (Tsukuba, 2000), JAERI-Conf
2001-002, pp. 357-376, in the proceedings https://doi.org/10.11484/jaeri-conf-2001-002, and open at
https://www.osti.gov/etdeweb/servlets/purl/20226979: the main reference, for eq. (1), the intensity,
eq. (3), the energy resolution, and eq. (4), the optimal ratio of the opening times of the two
choppers.

A. M. Gaspar, "Methods for analytically estimating the resolution and intensity of neutron
time-of-flight spectrometers. The case of the TOFTOF spectrometer", arXiv:0710.5319 (2007), a
technical report of 2005, https://arxiv.org/abs/0710.5319: eqs. (1), (14), (19) and (20), cited by
number in the steps above.

N. Violini, J. Voigt, S. Pasini and T. Brückel, "A method to compute the covariance matrix of
wavevector-energy transfer for neutron time-of-flight spectrometers", Nucl. Instrum. Methods Phys.
Res. A 736 (2014) 31-39, https://doi.org/10.1016/j.nima.2013.10.042: the covariance of Q and E for
direct-geometry spectrometers with pulse-shaping and monochromatic choppers, by the group that
proposed T-Rex. Not used here; the next step beyond the energy width alone.
"""

import warnings
from dataclasses import dataclass

import numpy as np
import scipp as sc
from scipy.constants import electron_volt, h, m_n
from scipy.optimize import OptimizeWarning, curve_fit
from scipy.signal import fftconvolve

from mcdgs.core import speed_from_wavelength
from mcdgs.instruments.trex.choppers import DISKS, chopper_timings
from mcdgs.instruments.trex.model import TRex

FWHM_PER_SIGMA = 2 * np.sqrt(2 * np.log(2))  # the full width at half maximum of a Gaussian
FIT_FRACTION = 0.2  # as the reduction: the Gaussian is fitted to the bins of at least this
GRID_STEPS = 400  # the narrowest of the three distributions spans this many steps of the grid
PAIRS = {"pulse_shaping": ("ps1", "ps2"), "monochromatic": ("mc1", "mc2")}


@dataclass(frozen=True)
class ChopperPair:
    """A counter-rotating pair of chopper disks, as the resolution sees it.

    Attributes
    ----------
    distance:
        From the moderator to the middle of the pair, along the beam, in m.
    spacing:
        From the first disk to the second, along the beam, in m.
    slit_arc:
        The width of the slit on the beam, at the rim, in m.
    rim_speed:
        The speed of the rim, in m/s, the same for both disks.
    beam_width:
        The width of the beam window of the disks, in m.
    """

    distance: float
    spacing: float
    slit_arc: float
    rim_speed: float
    beam_width: float


@dataclass(frozen=True)
class ResolutionGeometry:
    """What the elastic resolution depends on.

    Attributes
    ----------
    pulse_shaping, monochromatic:
        The two chopper pairs.
    sample_distance:
        From the middle of the monochromatic pair to the sample, ``L_MS``, in m.
    detector_distance:
        From the sample to the detector, ``L2``, in m.
    sample_radius:
        The radius of the cylinder of the sample, in m.
    two_theta_ranges:
        The scattering angles the detector sees, as (start, stop) in deg.
    lambda0:
        The wavelength the choppers are phased to, in Å.
    """

    pulse_shaping: ChopperPair
    monochromatic: ChopperPair
    sample_distance: float
    detector_distance: float
    sample_radius: float
    two_theta_ranges: tuple[tuple[float, float], ...]
    lambda0: float

    @property
    def coefficients(self) -> tuple[float, float]:
        """Return ``A`` and ``B`` of step 1, the weights of the two offsets in the flight time.

        ``A`` weighs the offset at the monochromatic pair, ``B`` that at the pulse-shaping pair.

        Returns
        -------
        :
            ``A = 1 + L_MD / L_PM`` and ``B = L_MD / L_PM``.
        """
        l_pm = self.monochromatic.distance - self.pulse_shaping.distance
        l_md = self.sample_distance + self.detector_distance
        return 1 + l_md / l_pm, l_md / l_pm


def read_pair(instrument: TRex, names: tuple[str, str]) -> ChopperPair:
    """Read a counter-rotating pair of the instrument, in its mode.

    Parameters
    ----------
    instrument:
        The instrument, with turning choppers.
    names:
        The two disks, in beam order, such as ``("ps1", "ps2")``.

    Returns
    -------
    :
        The pair, with the slit of the mode of ``instrument.params``.

    Raises
    ------
    ValueError
        If a disk is parked, if the disks turn at different speeds, or if their rims do not move
        across the beam in opposite directions.
    """
    params = instrument.params
    timings = chopper_timings(
        params.lambda0, params.rrm, instrument.source.frequency, params.mode, params.t0_offset
    )  # as the instrument computes them for McStas
    nu = {f"nu_{name}": timing.nu for name, timing in timings.items()}
    first, second = (DISKS[name] for name in names)
    # the rim moves across the beam one way or the other: by the sense of rotation, and by the
    # side of the disk the beam is on, beam_angle 0 or 180 deg
    motion = [
        np.sign(nu[f"nu_{name}"]) * np.cos(np.deg2rad(DISKS[name].beam_angle)) for name in names
    ]
    if 0 in motion:
        raise ValueError(f"{names}: a parked disk has no opening time")
    if abs(nu[f"nu_{names[0]}"]) != abs(nu[f"nu_{names[1]}"]) or first.radius != second.radius:
        raise ValueError(f"{names}: the disks do not turn at the same rim speed")
    if motion[0] == motion[1]:
        raise ValueError(
            f"{names}: the rims move across the beam the same way, not counter-rotating"
        )
    mode = params.mode
    _, slit_deg = first.slits[first.mode_slit[mode]]
    _, second_deg = second.slits[second.mode_slit[mode]]
    if slit_deg != second_deg:
        raise ValueError(f"{names}: the slits on the beam are {slit_deg} and {second_deg} deg")
    return ChopperPair(
        distance=(first.distance + second.distance) / 2,
        spacing=second.distance - first.distance,
        slit_arc=float(np.deg2rad(slit_deg) * first.radius),
        rim_speed=float(2 * np.pi * abs(nu[f"nu_{names[0]}"]) * first.radius),
        beam_width=float(instrument.choppers[names[0]].xwidth),
    )


def geometry_of(
    instrument: TRex,
    detector_distance: float,
    sample_radius: float,
    two_theta_ranges: tuple[tuple[float, float], ...],
) -> ResolutionGeometry:
    """Read the geometry of the resolution from the instrument, in its mode.

    Parameters
    ----------
    instrument:
        The instrument, with turning choppers and ``sample_pos``.
    detector_distance:
        From the sample to the detector, ``L2``, in m.
    sample_radius:
        The radius of the cylinder of the sample, in m.
    two_theta_ranges:
        The scattering angles the detector sees, as (start, stop) in deg.

    Returns
    -------
    :
        The geometry.
    """
    pairs = {name: read_pair(instrument, disks) for name, disks in PAIRS.items()}
    positions = instrument.positions()
    monochromatic = np.mean([positions[name] for name in PAIRS["monochromatic"]], axis=0)
    lambda0 = instrument.params.lambda0
    if lambda0 is None:
        raise ValueError("lambda0 is None: the white beam has no sub-pulses")
    return ResolutionGeometry(
        pulse_shaping=pairs["pulse_shaping"],
        monochromatic=pairs["monochromatic"],
        sample_distance=float(np.linalg.norm(np.array(positions["sample_pos"]) - monochromatic)),
        detector_distance=detector_distance,
        sample_radius=sample_radius,
        two_theta_ranges=two_theta_ranges,
        lambda0=lambda0,
    )


def pair_transmission(pair: ChopperPair, offset: float, times: np.ndarray) -> np.ndarray:
    """Compute how much of the beam passes a counter-rotating pair, at each time (step 2).

    The center of the slit of the first disk crosses the beam at ``u t``, that of the second at
    ``-u (t - offset)``. At ``x`` across the beam, of width ``w``, a neutron passes if it is inside
    both slits; the transmission is the length of the beam where it does, over ``w``.

    Parameters
    ----------
    pair:
        The pair.
    offset:
        The time from when the first slit is centered on the beam to when the second is, for the
        neutrons of the sub-pulse, in s.
    times:
        The times, in s.

    Returns
    -------
    :
        The fraction of the beam that passes, at each time.

    Examples
    --------
    Slits of 20 mm moving at 400 m/s over a beam of 30 mm are open as a triangle of 25 µs FWHM,
    ``a / (2 u)``, eq. (1) of Gaspar (2007), and pass at most 20 mm of the 30:

    >>> pair = ChopperPair(distance=0, spacing=0, slit_arc=0.02, rim_speed=400, beam_width=0.03)
    >>> opening = pair_transmission(pair, 0.0, np.array([0.0, 12.5e-6, 25e-6]))
    >>> [round(float(value), 6) for value in opening]
    [0.666667, 0.333333, 0.0]
    """
    u, a, w = pair.rim_speed, pair.slit_arc, pair.beam_width
    first = (u * times - a / 2, u * times + a / 2)
    second = (-u * (times - offset) - a / 2, -u * (times - offset) + a / 2)
    low = np.maximum.reduce([first[0], second[0], np.full_like(times, -w / 2)])
    high = np.minimum.reduce([first[1], second[1], np.full_like(times, w / 2)])
    return np.clip(high - low, 0, None) / w


def sample_path_density(
    radius: float, two_theta_ranges: tuple[tuple[float, float], ...], paths: np.ndarray
) -> np.ndarray:
    """Compute the distribution of the extra flight path from the sample (step 3).

    For one 2θ, the extra path of a uniform disk of ``radius`` is a semicircle of radius
    ``2 radius |sin θ|``; it is averaged over ``two_theta_ranges``, evenly per degree.

    Parameters
    ----------
    radius:
        The radius of the cylinder, in m.
    two_theta_ranges:
        The scattering angles, as (start, stop) in deg.
    paths:
        The extra paths, evenly spaced, in m.

    Returns
    -------
    :
        The density at each path, in 1/m.

    Examples
    --------
    It is a distribution: it sums to 1 over the paths.

    >>> paths = np.linspace(-0.03, 0.03, 6001)
    >>> density = sample_path_density(0.01, ((5, 144),), paths)
    >>> round(float(density.sum() * (paths[1] - paths[0])), 3)
    1.0
    """
    angles = np.concatenate([np.linspace(start, stop, 721) for start, stop in two_theta_ranges])
    weights = np.concatenate(
        [np.full(721, (stop - start) / 721) for start, stop in two_theta_ranges]
    )
    density = np.zeros_like(paths)
    for angle, weight in zip(angles, weights):
        reach = 2 * radius * abs(np.sin(np.deg2rad(angle) / 2))
        if reach == 0:
            continue
        inside = np.clip(1 - (paths / reach) ** 2, 0, None)
        density += weight * 2 * np.sqrt(inside) / (np.pi * reach)  # a semicircle of area 1
    return density / weights.sum()


def _gaussian(x: np.ndarray, height: float, center: float, sigma: float) -> np.ndarray:
    return height * np.exp(-((x - center) ** 2) / (2 * sigma**2))


def fwhm_of(energies: np.ndarray, density: np.ndarray, fit_fraction: float = FIT_FRACTION) -> float:
    """Fit a Gaussian to the core of a distribution, as the reduction, and return its FWHM (step 6).

    Parameters
    ----------
    energies:
        The energies, evenly spaced, in meV.
    density:
        The distribution at each energy.
    fit_fraction:
        The lowest value fitted, as a fraction of the highest.

    Returns
    -------
    :
        ``FWHM_PER_SIGMA`` times the σ of the Gaussian, in meV.

    Examples
    --------
    >>> x = np.linspace(-5, 5, 1001)
    >>> round(fwhm_of(x, np.exp(-(x**2) / 2)), 4)
    2.3548
    """
    core = density >= fit_fraction * density.max()
    x, y = energies[core], density[core]
    mean = np.average(x, weights=y)
    start = (y.max(), mean, np.sqrt(np.average((x - mean) ** 2, weights=y)))
    # weighted as the reduction weighs a bin of counts: by an error of the square root of it. The
    # covariance is not used: a line that is a Gaussian, as in the example, fits it exactly and
    # leaves it undefined, which scipy warns of
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", OptimizeWarning)
        best, _ = curve_fit(_gaussian, x, y, p0=start, sigma=np.sqrt(y))
    return float(FWHM_PER_SIGMA * abs(best[2]))


def opening_of(pair: ChopperPair, offset: float = 0.0) -> tuple[float, float]:
    """Compute the FWHM and the area of the opening of a pair, from :func:`pair_transmission`.

    Parameters
    ----------
    pair:
        The pair.
    offset:
        As :func:`pair_transmission`, in s.

    Returns
    -------
    :
        The full width at half maximum of the opening, in s, and its area, the fraction of the beam
        that passes integrated over time, in s.

    Examples
    --------
    Over a beam at least as wide as the slit, the triangle of eq. (1) of Gaspar (2007):

    >>> pair = ChopperPair(distance=0, spacing=0, slit_arc=0.02, rim_speed=400, beam_width=0.03)
    >>> fwhm, area = opening_of(pair)
    >>> round(fwhm * 1e6, 2), round(area * 1e6, 2)
    (25.0, 16.67)
    """
    reach = (pair.slit_arc + pair.beam_width) / (2 * pair.rim_speed) + abs(offset)
    times = np.linspace(-reach, reach, 200001)
    opening = pair_transmission(pair, offset, times)
    above = times[opening >= opening.max() / 2]
    return float(above[-1] - above[0]), float(opening.sum() * (times[1] - times[0]))


def slit_opening(pair: ChopperPair) -> float:
    """Return the FWHM of the opening of a pair from its slit alone, ``a / (2 u)``.

    Eq. (1) of Gaspar (2007), ``Δt = b / (2π f D)``: the triangle of a beam at least as wide as
    the slit, whatever the beam.

    Parameters
    ----------
    pair:
        The pair.

    Returns
    -------
    :
        The FWHM, in s.
    """
    return pair.slit_arc / (2 * pair.rim_speed)


def sample_path_spread(radius: float, two_theta_ranges: tuple[tuple[float, float], ...]) -> float:
    """Return the spread of the flight path from the sample of eq. (20) of Gaspar (2007).

    On the elastic line, eq. (20) is ``v Δt_s = sqrt(c² sin² 2θ + d² (1 - cos 2θ)²)``, for a
    rectangle of half-width ``c`` across the beam and half-thickness ``d`` along it, which Gaspar
    calls ``a``, here kept for the arc of a slit. For the cylinder, both are its ``radius``. It is averaged over ``two_theta_ranges`` as a root mean
    square, evenly per degree.

    Parameters
    ----------
    radius:
        The radius of the cylinder, in m.
    two_theta_ranges:
        The scattering angles, as (start, stop) in deg.

    Returns
    -------
    :
        ``ΔL``, the spread of the flight path, in m.

    Examples
    --------
    At 90°, a square of half-width 1 cm spreads the path by the root of 2 cm²:

    >>> round(sample_path_spread(0.01, ((89.999, 90.001),)) * 100, 4)
    1.4142
    """
    angles = np.deg2rad(
        np.concatenate([np.linspace(start, stop, 721) for start, stop in two_theta_ranges])
    )
    weights = np.concatenate(
        [np.full(721, (stop - start) / 721) for start, stop in two_theta_ranges]
    )
    spread = radius**2 * (np.sin(angles) ** 2 + (1 - np.cos(angles)) ** 2)
    return float(np.sqrt(np.average(spread, weights=weights)))


def literature_resolution(
    geometry: ResolutionGeometry, incident_energies: sc.Variable, openings: str = "slit"
) -> sc.Dataset:
    """Compute the FWHM of the elastic line from eq. (3) of Lechner (2000), on the elastic line.

    ``Δ(hν) = (2 Ei v0 / (L_PM L2)) sqrt((τ1 L_MD)² + (τ2 (L_PM + L_MD))² + (L_PM ΔL / v0)²)``:
    eq. (3) at ``λ = λ0``, with ``2 Ei v0`` for his 647.2 / λ³ in µeV, Å, m and µs, ``τ1`` and
    ``τ2`` the FWHM of the openings of the pulse-shaping and monochromatic pairs, and ``ΔL`` the
    spread of the flight path, here from the sample only, :func:`sample_path_spread`. The terms add
    in quadrature, as for Gaussians.

    Parameters
    ----------
    geometry:
        The geometry, as :func:`geometry_of` reads it.
    incident_energies:
        The incident energies, along one dimension.
    openings:
        ``"slit"`` for the opening times of the slits, :func:`slit_opening`; ``"beam"`` for
        the exact FWHM of the openings for the width of the beam, :func:`opening_of`.

    Returns
    -------
    :
        Along the dimension of ``incident_energies``, in meV: ``fwhm``, and ``pulse_shaping``,
        ``monochromatic`` and ``sample``, the three terms; with ``incident_energy`` as coordinate.

    Raises
    ------
    ValueError
        If ``openings`` is neither ``"slit"`` nor ``"beam"``.
    """
    if openings not in ("slit", "beam"):
        raise ValueError(f"openings {openings!r} is not 'slit' or 'beam'")
    dim = incident_energies.dim
    energies = incident_energies.to(unit="meV").values
    l_pm = geometry.monochromatic.distance - geometry.pulse_shaping.distance
    l_md = geometry.sample_distance + geometry.detector_distance
    spread = sample_path_spread(geometry.sample_radius, geometry.two_theta_ranges)
    v_lambda0 = speed_from_wavelength(geometry.lambda0)
    widths: dict[str, list[float]] = {
        name: [] for name in ("fwhm", "pulse_shaping", "monochromatic", "sample")
    }
    for ei in energies:
        v0 = _speed(float(ei))
        taus = {}
        for name in PAIRS:
            pair: ChopperPair = getattr(geometry, name)
            if openings == "slit":
                taus[name] = slit_opening(pair)
            else:
                taus[name], _ = opening_of(pair, pair.spacing * (1 / v0 - 1 / v_lambda0))
        scale = 2 * ei * v0 / (l_pm * geometry.detector_distance)
        terms = {
            "pulse_shaping": scale * taus["pulse_shaping"] * l_md,
            "monochromatic": scale * taus["monochromatic"] * (l_pm + l_md),
            "sample": scale * l_pm * spread / v0,
        }
        for name, value in terms.items():
            widths[name].append(float(value))
        widths["fwhm"].append(float(np.sqrt(sum(value**2 for value in terms.values()))))
    return sc.Dataset(
        {name: sc.array(dims=[dim], values=values, unit="meV") for name, values in widths.items()},
        coords={"incident_energy": incident_energies.to(unit="meV")},
    )


def optimal_ratio(geometry: ResolutionGeometry) -> float:
    """Return the optimal ratio of the opening times of the two pairs, eq. (4) of Lechner (2000).

    On the elastic line, ``ρ_opt = τ1 / τ2 = (L_PM + L_MD) / L_MD``, ``A / B`` of step 1: the
    ratio at which the two chopper terms are equal, which gives the most intensity for a resolution.

    Parameters
    ----------
    geometry:
        The geometry.

    Returns
    -------
    :
        ``ρ_opt``.
    """
    a_coefficient, b_coefficient = geometry.coefficients
    return a_coefficient / b_coefficient


def energy_of(wavelength: float) -> float:
    """Return the energy of a neutron of ``wavelength`` Å, in meV.

    Examples
    --------
    >>> round(energy_of(2.0), 2)
    20.45
    """
    return float((h / (wavelength * 1e-10)) ** 2 / (2 * m_n) / electron_volt * 1e3)


def _speed(energy: float) -> float:
    """Return the speed of a neutron of ``energy`` meV, in m/s."""
    return float(np.sqrt(2 * energy * 1e-3 * electron_volt / m_n))


def line_shapes(
    geometry: ResolutionGeometry, ei: float
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Compute the distribution of δE of each term, and of the three together, at one Ei.

    Steps 1 to 5: each term on one grid of δE, and their convolution.

    Parameters
    ----------
    geometry:
        The geometry, as :func:`geometry_of` reads it.
    ei:
        The incident energy, in meV.

    Returns
    -------
    :
        By term, ``pulse_shaping``, ``monochromatic``, ``sample`` and ``total``: the δE of the grid,
        in meV, and the distribution there, summing to 1 for the three terms.

    Examples
    --------
    The three terms together are wider than each alone:

    >>> from dataclasses import replace
    >>> pair = ChopperPair(distance=0, spacing=0, slit_arc=0.1, rim_speed=400, beam_width=0.05)
    >>> geometry = ResolutionGeometry(
    ...     pulse_shaping=pair,
    ...     monochromatic=replace(pair, distance=50, slit_arc=0.02),
    ...     sample_distance=2,
    ...     detector_distance=3,
    ...     sample_radius=0.01,
    ...     two_theta_ranges=((5, 140),),
    ...     lambda0=2,
    ... )
    >>> widths = {name: fwhm_of(*line) for name, line in line_shapes(geometry, 20).items()}
    >>> widths["total"] > max(widths[name] for name in ("pulse_shaping", "monochromatic", "sample"))
    True
    """
    v0 = _speed(ei)
    k = 2 * ei * v0 / geometry.detector_distance  # meV per s of flight time, step 1
    a_coefficient, b_coefficient = geometry.coefficients
    v_lambda0 = speed_from_wavelength(geometry.lambda0)
    # the energy of each term per its offset: δE = scale × offset
    scales = {
        "pulse_shaping": -k * b_coefficient,
        "monochromatic": k * a_coefficient,
        "sample": 2 * ei / geometry.detector_distance,
    }
    reaches = {}  # the largest offset of each term: in s for the pairs, in m for the sample
    offsets = {}
    for name in PAIRS:
        pair: ChopperPair = getattr(geometry, name)
        offsets[name] = pair.spacing * (1 / v0 - 1 / v_lambda0)
        reaches[name] = (pair.slit_arc + pair.beam_width) / (2 * pair.rim_speed) + abs(
            offsets[name]
        )
    reaches["sample"] = 2 * geometry.sample_radius
    spans = {name: abs(scales[name]) * reach for name, reach in reaches.items()}
    step = min(spans.values()) / GRID_STEPS
    half = sum(spans.values())
    energies = np.arange(-half, half + step / 2, step)
    terms = {}
    for name in PAIRS:
        times = energies / scales[name]
        # the opening is centered on its mean: its width does not depend on where it is
        opening = pair_transmission(getattr(geometry, name), offsets[name], times)
        terms[name] = opening / opening.sum()
    path_density = sample_path_density(
        geometry.sample_radius, geometry.two_theta_ranges, energies / scales["sample"]
    )
    terms["sample"] = path_density / path_density.sum()
    total = terms["pulse_shaping"]
    for name in ("monochromatic", "sample"):
        total = fftconvolve(total, terms[name], mode="same")
    terms["total"] = np.clip(total, 0, None)
    return {name: (energies, density) for name, density in terms.items()}


def elastic_resolution(
    geometry: ResolutionGeometry, incident_energies: sc.Variable, fit_fraction: float = FIT_FRACTION
) -> sc.Dataset:
    """Compute the FWHM of the elastic line at each incident energy, and that of each term.

    Parameters
    ----------
    geometry:
        The geometry, as :func:`geometry_of` reads it.
    incident_energies:
        The incident energies, along one dimension.
    fit_fraction:
        As :func:`fwhm_of`.

    Returns
    -------
    :
        Along the dimension of ``incident_energies``, in meV: ``fwhm``, of the three terms
        together, and ``pulse_shaping``, ``monochromatic`` and ``sample``, the FWHM of each term
        alone; with ``incident_energy`` as coordinate.
    """
    dim = incident_energies.dim
    energies = incident_energies.to(unit="meV").values
    widths: dict[str, list[float]] = {
        name: [] for name in ("fwhm", "pulse_shaping", "monochromatic", "sample")
    }
    for ei in energies:
        lines = line_shapes(geometry, float(ei))
        widths["fwhm"].append(fwhm_of(*lines["total"], fit_fraction))
        for name in ("pulse_shaping", "monochromatic", "sample"):
            widths[name].append(fwhm_of(*lines[name], fit_fraction))
    return sc.Dataset(
        {name: sc.array(dims=[dim], values=values, unit="meV") for name, values in widths.items()},
        coords={"incident_energy": incident_energies.to(unit="meV")},
    )
