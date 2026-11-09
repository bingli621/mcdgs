"""Positions in the beam frame of the sample.

The momentum transfer of a neutron is ``ki - kf``, along the incident beam and across it, so
its components need a frame tied to the beam. :func:`beam_frame` gives positions in it: the
origin is the sample, z points along the incident beam, y up and x to the left of the beam, as
in McStas. A guide that bends the beam tilts this frame from the lab frame.
"""

import scipp as sc

__all__ = ["beam_frame"]


def beam_frame(
    position: sc.Variable, sample_position: sc.Variable, upstream_position: sc.Variable
) -> sc.Variable:
    """Return positions in the beam frame of the sample.

    Its origin is the sample, z points along the beam, from ``upstream_position`` to the sample,
    y up, as the y of the lab, and x to the left of the beam, ``y × z``, as in McStas.

    Parameters
    ----------
    position:
        The positions, as vectors in the lab frame.
    sample_position:
        Where the sample is.
    upstream_position:
        A point on the beam before the sample, such as a monitor.

    Returns
    -------
    :
        The positions in the beam frame, in the unit of ``position``.

    Examples
    --------
    >>> local = beam_frame(
    ...     sc.vector([2.0, 0.0, 1.0], unit="m"),
    ...     sc.vector([1.0, 0.0, 1.0], unit="m"),
    ...     sc.vector([0.0, 0.0, 0.0], unit="m"),
    ... )
    >>> [round(value, 4) for value in local.value.tolist()]  # 45° to the left of a beam at 45°
    [0.7071, 0.0, 0.7071]
    """
    beam = sample_position - upstream_position.to(unit=sample_position.unit)
    z = beam / sc.norm(beam)
    y = sc.vector([0.0, 1.0, 0.0])
    x = sc.cross(y, z)
    x = x / sc.norm(x)
    y = sc.cross(z, x)  # exactly across the beam, should it not be horizontal
    relative = position - sample_position.to(unit=position.unit)
    return sc.spatial.as_vectors(sc.dot(relative, x), sc.dot(relative, y), sc.dot(relative, z))
