"""Powders that scatter inelastically: a ferromagnet, from one magnon branch.

:func:`magnon_powder` returns an :class:`.Isotropic_Sqw` that reads the S(|Q|, E) of a file, which
:func:`write_magnon_sqw` writes: the magnon of a simple-cubic ferromagnet, averaged over the
directions of Q, as a powder has it.
"""

from pathlib import Path
from typing import Any

import numpy as np
from scipy.ndimage import gaussian_filter1d

from mcdgs.core.components import Isotropic_Sqw

__all__ = [
    "magnon_band",
    "magnon_energies",
    "magnon_powder",
    "sphere_directions",
    "write_magnon_sqw",
]

# the default magnon: a gap of MAGNON_GAP at Q = 0, up to MAGNON_TOP at the corner of the zone
MAGNON_A = 4.0  # Å, the lattice constant of the simple-cubic cell
MAGNON_GAP = 5.0  # meV
MAGNON_TOP = 20.0  # meV
MAGNON_JS = (MAGNON_TOP - MAGNON_GAP) / 24  # meV: the band is 24 JS wide
MAGNON_WIDTH = 0.03  # meV, the σ of the line: far below the energy resolution of T-Rex
# the S(|Q|, E) table: its axes, evenly spaced, and the directions it averages over
SQW_Q_MAX = 10.0  # 1/Å, past the largest |Q| of T-Rex for E of MAGNON_GAP or more
SQW_Q_STEP = 0.02  # 1/Å
SQW_E_STEP = 0.02  # meV
SQW_DIRECTIONS = 20_000
MAGNON_SIGMA_COH = 5.0  # barns per cell: only the scale of the intensity


def sphere_directions(n: int = SQW_DIRECTIONS) -> np.ndarray:
    """Return ``n`` unit vectors spread evenly over the sphere, on a Fibonacci lattice.

    Unlike random directions, they are the same at every call, so the table is too.

    Examples
    --------
    >>> directions = sphere_directions(1000)
    >>> directions.shape, bool(np.allclose(np.linalg.norm(directions, axis=1), 1))
    ((1000, 3), True)
    >>> bool(np.allclose(directions.mean(axis=0), 0, atol=1e-3))  # no direction is favored
    True
    """
    k = np.arange(n) + 0.5
    z = 1 - 2 * k / n
    phi = np.pi * (1 + np.sqrt(5)) * k  # the golden angle, times k
    r = np.sqrt(1 - z**2)
    return np.column_stack([r * np.cos(phi), r * np.sin(phi), z])


def magnon_energies(
    q: np.ndarray,
    a: float = MAGNON_A,
    js: float = MAGNON_JS,
    gap: float = MAGNON_GAP,
) -> np.ndarray:
    """Return the energy of the magnon of a simple-cubic ferromagnet at each wavevector.

    With an exchange J between nearest neighbors, spins S and an anisotropy gap, it is
    ``gap + 4 J S (3 - cos(qx a) - cos(qy a) - cos(qz a))``: ``gap`` at the center of the zone,
    ``gap + 24 J S`` at its corner, (π/a)(1, 1, 1).

    Parameters
    ----------
    q:
        Wavevectors, in 1/Å, along a last axis of three.
    a:
        Lattice constant, in Å.
    js:
        Exchange times spin, J S, in meV.
    gap:
        Energy at the center of the zone, in meV.

    Returns
    -------
    :
        The energies, in meV, with the shape of ``q`` without its last axis.

    Examples
    --------
    >>> corner = np.full(3, np.pi / MAGNON_A)
    >>> float(magnon_energies(np.zeros(3))), round(float(magnon_energies(corner)), 9)
    (5.0, 20.0)
    """
    return gap + 4 * js * (3 - np.cos(np.asarray(q) * a).sum(axis=-1))


def magnon_band(
    q: np.ndarray,
    a: float = MAGNON_A,
    js: float = MAGNON_JS,
    gap: float = MAGNON_GAP,
    directions: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the lowest and the highest energy of the magnon at each |Q|, over its directions.

    A powder scatters at a |Q| from every wavevector of that length, so at a |Q| the magnon of a
    powder spans this band. Its edges, sharp in the model, are where the resolution shows.

    Parameters
    ----------
    q:
        The values of |Q|, in 1/Å.
    a, js, gap:
        As :func:`magnon_energies`.
    directions:
        Unit vectors to average over. If ``None``, those of :func:`sphere_directions`.

    Returns
    -------
    :
        The lowest and the highest energy at each |Q|, in meV.

    Examples
    --------
    >>> low, high = magnon_band(np.array([0.0]))
    >>> float(low[0]), float(high[0])
    (5.0, 5.0)
    """
    directions = sphere_directions() if directions is None else directions
    low, high = np.empty(len(q)), np.empty(len(q))
    for i, length in enumerate(q):
        energies = magnon_energies(length * directions, a, js, gap)
        low[i], high[i] = energies.min(), energies.max()
    return low, high


def write_magnon_sqw(
    path: str | Path,
    a: float = MAGNON_A,
    js: float = MAGNON_JS,
    gap: float = MAGNON_GAP,
    width: float = MAGNON_WIDTH,
    q_max: float = SQW_Q_MAX,
    q_step: float = SQW_Q_STEP,
    e_step: float = SQW_E_STEP,
    directions: int = SQW_DIRECTIONS,
) -> Path:
    """Write the S(|Q|, E) of the magnon of a powder, for :class:`.Isotropic_Sqw`.

    At each |Q| of the table, the magnon of :func:`magnon_energies` is averaged over the
    ``directions`` of :func:`sphere_directions`: the energies they give are histogrammed and
    smoothed by a Gaussian of σ ``width``, the line of the magnon. Each row is a density in E:
    the sum of each side times ``e_step`` is 1, one magnon, the same at every |Q|. The magnetic
    form factor and the polarization factor are left out, so the intensity does not fall with |Q|.

    E is the energy the neutron loses. The magnon is created at positive E and annihilated at
    negative E, as often: the table is symmetric in E, from a little below the bottom of the band
    on the gain side to a little above its top on the loss side, a Bose factor of 1, the limit of
    a high temperature. It holds both sides itself, rather than leaving :class:`.Isotropic_Sqw` to
    mirror one, and ``classical`` 0 in its header: with ``T=-1``, no detailed balance, the powder
    scatters it as it is.

    Parameters
    ----------
    path:
        The file to write; its folder is created.
    a, js, gap:
        As :func:`magnon_energies`.
    width:
        σ of the line, in meV.
    q_max, q_step:
        The values of |Q|, from 0 to ``q_max``, in 1/Å.
    e_step:
        The spacing of E, in meV.
    directions:
        The number of directions averaged over.

    Returns
    -------
    :
        The path of the file.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    q = np.arange(0.0, q_max + q_step / 2, q_step)
    e_max = gap + 24 * js + 6 * width
    loss = np.arange(0.0, e_max + e_step / 2, e_step)  # E from 0 up: the magnon created
    edges = np.append(loss - e_step / 2, loss[-1] + e_step / 2)
    unit = sphere_directions(directions)
    created = np.empty((q.size, loss.size))
    for i, length in enumerate(q):
        created[i] = np.histogram(magnon_energies(length * unit, a, js, gap), bins=edges)[0]
    created = gaussian_filter1d(created / (directions * e_step), width / e_step, axis=1)
    # the gain side, the magnon annihilated: the loss side mirrored, a Bose factor of 1
    energy = np.concatenate([-loss[:0:-1], loss])
    table = np.concatenate([created[:, :0:-1], created], axis=1)
    with path.open("w") as f:
        f.write(
            "# S(|Q|, E) of the magnon of a simple-cubic ferromagnet, averaged over the "
            "directions of Q\n"
            f"# a = {a} Å, J S = {js} meV, gap = {gap} meV, line σ = {width} meV, "
            f"{directions} directions\n"
            "# classical 0\n"
            "# q axis values [1/Å]\n"
        )
        np.savetxt(f, q[None], fmt="%.4f")
        f.write("# w axis values [meV], the energy the neutron loses, negative when it gains\n")
        np.savetxt(f, energy[None], fmt="%.4f")
        f.write("# sqw values [1/meV], one row per q value\n")
        np.savetxt(f, table, fmt="%.4e")
    return path


def magnon_powder(sqw_file: str | Path, **params: Any) -> Isotropic_Sqw:
    """Create a powder that scatters from the magnon of the S(|Q|, E) of ``sqw_file``.

    A solid cylinder 20 mm across and 10 mm high, the size of :func:`.vanadium`, so the sample
    adds the same to the resolution. One simple-cubic cell of :data:`MAGNON_A` is the
    scattering unit, with a coherent cross-section of :data:`MAGNON_SIGMA_COH`, which sets only
    the scale of the intensity; it neither absorbs nor scatters incoherently. S(|Q|, E) is taken
    as it is: no detailed balance (``T=-1``, ``classical=0``) and no scaling (``norm=1``).

    No focusing is set: it depends on the instrument. Give the height of the detector seen from
    the sample with ``d_phi`` and ``order=1``.

    Parameters
    ----------
    sqw_file:
        The S(|Q|, E), such as :func:`write_magnon_sqw` writes. A relative path is one from the
        folder McStas runs in, so give an absolute one.
    **params:
        Further ``Isotropic_Sqw`` parameters, or changes to the defaults (e.g. ``name``,
        ``at``, ``relative``, ``order``, ``d_phi``, ``p_interact``).

    Returns
    -------
    :
        The sample.
    """
    return Isotropic_Sqw(
        **{
            "Sqw_coh": str(sqw_file),
            "radius": 0.01,
            "yheight": 0.01,
            "rho": 1 / MAGNON_A**3,  # 1/Å³, cells
            "sigma_coh": MAGNON_SIGMA_COH,  # barns
            "sigma_inc": -1,  # none: an empty Sqw_inc would scatter like vanadium
            "sigma_abs": 0,
            "T": -1,
            "classical": 0,
            "norm": 1,
            **params,
        }
    )
