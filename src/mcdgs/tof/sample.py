"""The tof sample of a McStas sample that changes the energy of the neutrons.

tof has no scattering: its :class:`tof.InelasticSample` gives every neutron that reaches it a new
energy, and sends it on along the same flight path. So :func:`inelastic_sample` keeps only the
energy transfers of the McStas sample, as two bands, the neutron losing energy and gaining it,
and draws each neutron's transfer uniformly over them. The time-distance diagram then shows every
arrival time at the detector that the sample can give, which is what a frame overlap needs, but
not the intensities.

``w`` is the energy the neutron loses, so its final energy is ``E_f = E_i - w``. A sample
scatters over a band of ``|w|``, from ``w_min`` to ``w_max`` (:func:`energy_transfer_band`), and
tof draws ``w`` from ``-w_max`` to ``-w_min`` and from ``w_min`` to ``w_max``, as often, whichever
sides the McStas sample has. A neutron that would end with a negative energy is removed by tof.
The band of each McStas sample that changes the energy:

- ``Isotropic_Sqw``: where its S(q, w) file scatters. The file has three blocks of numbers: the
  values of q, those of w, and S(q, w), one row per q. McStas leaves out the values of S(q, w)
  below ``threshold`` (``Isotropic_Sqw.comp``, ``Sqw_readfile``), so the band runs from the lowest
  to the highest ``|w|`` where any S(q, w) is at least ``threshold`` and above 0. The sample takes
  the widest band of its coherent and its incoherent file. A file of one block is a powder file
  of Bragg lines, and a part with no S(q, w) above ``threshold`` adds nothing. A part left out
  (``sigma_coh`` or ``sigma_inc`` -1) adds nothing either; the incoherent part with no file is
  elastic, a transfer of 0.
- ``Phonon_simple``: its branch, which creates and annihilates phonons, ``sqrt(gap**2 + (12 -
  J(q)) c**2 / a**2)`` (``Phonon_simple.comp``, ``omega_q``): from ``gap`` at q = 0 to
  ``sqrt(gap**2 + (4 c / a)**2)`` at the zone boundary along a cubic axis.
- ``Incoherent``: with a ``deltaE``, ``w`` is uniform from ``Etrans - |deltaE|`` to
  ``Etrans + |deltaE|`` (``Incoherent.comp``, TRACE). Without, the fraction ``f_QE`` of the
  scattering is quasielastic, a Lorentzian of half width ``gamma``, which has no end: its band is
  cut at ``QUASIELASTIC_GAMMAS`` times ``gamma``, which holds 94 % of it. The rest is elastic.

A sample that scatters only elastically, such as vanadium, or a component that is not a sample,
has no band and no tof sample: tof takes it as not there.
"""

from collections.abc import Callable
from pathlib import Path

import numpy as np
import scipp as sc
import tof

from mcdgs.core.components import Component, Incoherent, Isotropic_Sqw, Phonon_simple

SQW_BLOCKS = 3  # the blocks of an S(q, w) file: q, w and S(q, w)
POWDER_BLOCKS = 1  # the one block of a powder file, its Bragg lines
MISSING = -1  # a cross-section of -1 leaves its part out
QUASIELASTIC_GAMMAS = 10  # the band of a Lorentzian, in half widths: 2 / π atan(10), 94 %


def read_blocks(path: str | Path) -> list[np.ndarray]:
    r"""Return the blocks of numbers of a McStas table file.

    A block is a run of lines of numbers; any other line, a comment, a blank or text, ends it,
    as McStas's ``Table_Read_Array`` reads them.

    Parameters
    ----------
    path:
        The file.

    Returns
    -------
    :
        The numbers of each block, flat, in the order of the file.

    Examples
    --------
    >>> import tempfile
    >>> folder = Path(tempfile.mkdtemp())
    >>> _ = (folder / "s.sqw").write_text("# q\n0 1\n# w\n-2 0 2\n# sqw\n1 1 1\n1 1 1\n")
    >>> [block.tolist() for block in read_blocks(folder / "s.sqw")]
    [[0.0, 1.0], [-2.0, 0.0, 2.0], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]]
    """
    blocks: list[np.ndarray] = []
    block: list[float] = []
    with Path(path).open() as f:
        for line in f:
            try:
                values = [float(word) for word in line.split()]
            except ValueError:
                values = []
            if values:
                block.extend(values)
            elif block:
                blocks.append(np.array(block))
                block = []
    if block:
        blocks.append(np.array(block))
    return blocks


def sqw_energy_transfers(path: str | Path, threshold: float) -> tuple[float, float] | None:
    """Return the lowest and the highest ``|w|`` at which an ``Isotropic_Sqw`` file scatters.

    Parameters
    ----------
    path:
        An S(q, w) file, of three blocks, or a powder file, of one.
    threshold:
        The ``threshold`` of the sample: McStas leaves out the S(q, w) below it.

    Returns
    -------
    :
        ``(w_min, w_max)``, in meV, the lowest and the highest ``|w|`` where any S(q, w) is at
        least ``threshold`` and above 0; ``(0, 0)`` for a powder file, which is elastic; ``None``
        if no S(q, w) is.

    Raises
    ------
    ValueError
        If the file has neither three blocks nor one, or its S(q, w) is not one value for each
        q and w.
    """
    blocks = read_blocks(path)
    if len(blocks) == POWDER_BLOCKS:
        return 0.0, 0.0
    if len(blocks) != SQW_BLOCKS:
        raise ValueError(
            f"{path}: {len(blocks)} blocks of numbers, not {SQW_BLOCKS} (q, w, S(q, w)) or "
            f"{POWDER_BLOCKS} (a powder file)"
        )
    q, w, sqw = blocks
    if sqw.size != q.size * w.size:
        raise ValueError(f"{path}: {sqw.size} values of S(q, w), not {q.size} q by {w.size} w")
    sqw = sqw.reshape(q.size, w.size)
    scatters = ((sqw >= threshold) & (sqw > 0)).any(axis=0)
    if not scatters.any():
        return None
    transfers = np.abs(w[scatters])
    return float(transfers.min()), float(transfers.max())


def sqw_band(sample: Isotropic_Sqw) -> tuple[float, float]:
    """Return the band of ``|w|`` over which an ``Isotropic_Sqw`` scatters.

    Parameters
    ----------
    sample:
        The sample. Its files are read from where they are named; a relative path is from the
        folder Python runs in.

    Returns
    -------
    :
        ``(w_min, w_max)``, in meV, the widest of its coherent and its incoherent parts.
        ``(0, 0)`` if it scatters elastically only, or not at all.
    """
    bands: list[tuple[float, float]] = []
    parts = ((sample.Sqw_coh, sample.sigma_coh, False), (sample.Sqw_inc, sample.sigma_inc, True))
    for file, sigma, incoherent in parts:
        if sigma == MISSING:
            continue
        if file:  # McStas takes "" as no file, as None
            band = sqw_energy_transfers(file, sample.threshold)
            if band is not None:
                bands.append(band)
        elif incoherent:  # no file: elastic, the same in every direction, as vanadium
            bands.append((0.0, 0.0))
    if not bands:
        return 0.0, 0.0
    return min(low for low, _ in bands), max(high for _, high in bands)


def phonon_band(sample: Phonon_simple) -> tuple[float, float]:
    """Return the band of ``|w|`` of the branch of a ``Phonon_simple``.

    Parameters
    ----------
    sample:
        The sample.

    Returns
    -------
    :
        ``(w_min, w_max)``, in meV: the branch at q = 0, ``gap``, and at the zone boundary along
        a cubic axis, ``sqrt(gap**2 + (4 c / a)**2)``. McStas takes a ``gap`` of 0 or less as none.

    Examples
    --------
    >>> sample = Phonon_simple(
    ...     radius=0.005, yheight=0.01, sigma_abs=0, sigma_inc=0, a=4, b=1, M=1, c=3, DW=1, T=1
    ... )
    >>> phonon_band(sample)
    (0.0, 3.0)
    >>> from dataclasses import replace
    >>> phonon_band(replace(sample, gap=4.0))
    (4.0, 5.0)
    """
    gap = float(max(sample.gap, 0))
    return gap, float(np.hypot(gap, sample.branch_maximum))


def incoherent_band(sample: Incoherent) -> tuple[float, float]:
    """Return the band of ``|w|`` over which an ``Incoherent`` scatters.

    Parameters
    ----------
    sample:
        The sample.

    Returns
    -------
    :
        ``(w_min, w_max)``, in meV: around ``Etrans`` with a ``deltaE``; else from 0 to
        ``QUASIELASTIC_GAMMAS`` times ``gamma`` if some of it is quasielastic; else ``(0, 0)``,
        elastic.

    Examples
    --------
    >>> incoherent_band(Incoherent(Etrans=10, deltaE=2))
    (8.0, 12.0)
    >>> incoherent_band(Incoherent(Etrans=1, deltaE=2))  # across the elastic line
    (0.0, 3.0)
    >>> incoherent_band(Incoherent(f_QE=0.5, gamma=0.1))
    (0.0, 1.0)
    >>> incoherent_band(Incoherent())  # vanadium
    (0.0, 0.0)
    """
    if sample.deltaE != 0:
        low, high = sample.Etrans - abs(sample.deltaE), sample.Etrans + abs(sample.deltaE)
        w_min = 0 if low <= 0 <= high else min(abs(low), abs(high))
        return float(w_min), float(max(abs(low), abs(high)))
    if sample.f_QE > 0 and sample.gamma != 0:
        return 0.0, float(QUASIELASTIC_GAMMAS * abs(sample.gamma))
    return 0.0, 0.0


def energy_transfer_band(component: Component) -> tuple[float, float] | None:
    """Return the band of ``|w|`` over which a McStas sample changes the energy, on either side.

    Parameters
    ----------
    component:
        Any component.

    Returns
    -------
    :
        ``(w_min, w_max)``, in meV: the sample scatters from ``-w_max`` to ``-w_min``, the
        neutron gaining energy, and from ``w_min`` to ``w_max``, losing it, through
        :func:`sqw_band`, :func:`phonon_band` or :func:`incoherent_band`. ``None`` for a sample
        that scatters elastically only, and for a component that is not one of them.
    """
    match component:
        case Isotropic_Sqw():
            band = sqw_band(component)
        case Phonon_simple():
            band = phonon_band(component)
        case Incoherent():
            band = incoherent_band(component)
        case _:
            return None
    return band if band[1] > 0 else None


def two_band_final_energy(
    w_min: float, w_max: float, rng: np.random.Generator
) -> Callable[[sc.Variable], sc.Variable]:
    """Return the final energy of tof neutrons, after a transfer drawn over two bands.

    Each neutron loses an energy ``w`` drawn uniformly from ``-w_max`` to ``-w_min`` and from
    ``w_min`` to ``w_max``: a ``|w|`` uniform from ``w_min`` to ``w_max``, of either sign as
    often.

    Parameters
    ----------
    w_min, w_max:
        The band of ``|w|``, in meV.
    rng:
        The random numbers, drawn anew at each call.

    Returns
    -------
    :
        The final energy ``E_f = E_i - w`` of each incident energy ``E_i``, in meV.

    Examples
    --------
    >>> final = two_band_final_energy(1.0, 2.0, np.random.default_rng(1))
    >>> ef = final(sc.full(dims=["event"], shape=[1000], value=10.0, unit="meV"))
    >>> transfer = np.abs(10 - ef.values)
    >>> bool(((transfer >= 1) & (transfer <= 2)).all()), bool((ef.values > 10).any())
    (True, True)
    """

    def final_energy(incident: sc.Variable) -> sc.Variable:
        incident = incident.to(unit="meV", dtype=float)
        size = incident.shape
        transfer = rng.uniform(w_min, w_max, size=size) * rng.choice([-1.0, 1.0], size=size)
        return incident - sc.array(dims=incident.dims, values=transfer, unit="meV")

    return final_energy


def inelastic_sample(
    name: str,
    component: Component,
    distance: float,
    seed: int | None = None,
) -> tof.InelasticSample | None:  # type: ignore
    """Return the tof sample of a McStas sample that changes the energy, over its two bands.

    Parameters
    ----------
    name:
        The name of the sample.
    component:
        The McStas component.
    distance:
        Its distance from the source, in m.
    seed:
        Seed of the random numbers of the transfers. They are a stream of their own, apart from
        those of the source, which take the same seed. If ``None``, every build differs.

    Returns
    -------
    :
        The tof sample, of :func:`two_band_final_energy` over :func:`energy_transfer_band`;
        ``None`` if the component does not change the energy.
    """
    band = energy_transfer_band(component)
    if band is None:
        return None
    w_min, w_max = band
    (sample_seed,) = np.random.SeedSequence(seed).spawn(1)  # not the source's own stream
    return tof.InelasticSample(
        distance=sc.scalar(distance, unit="m"),
        name=name,
        func=two_band_final_energy(w_min, w_max, np.random.default_rng(sample_seed)),
    )  # type: ignore
