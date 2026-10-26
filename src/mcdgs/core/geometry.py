"""Positions and orientations of components in the lab frame, as McStas computes them.

A component sits at ``at`` (m) in the frame of the component named ``relative`` (``None`` means
the lab frame). It is rotated by ``rotated`` (deg, about x, then y, then z) relative to that
frame, as McStas AT and ROTATED with RELATIVE.

A value may be a McStas expression of instrument parameters such as ``"-chi"`` (see
:func:`.utils.evaluate`). A rotation is the matrix from the lab frame to the component's own
frame, as ``rot_set_rotation`` in McStas's mccode-r.c.
"""

from collections.abc import Mapping
from typing import TYPE_CHECKING

import numpy as np

from .utils import evaluate

if TYPE_CHECKING:
    from .components import Component


def rotation(rx: float, ry: float, rz: float) -> np.ndarray:
    """Compute the McStas rotation matrix of ROTATED (rx, ry, rz).

    Parameters
    ----------
    rx, ry, rz:
        Rotation angles about x, y and z, in deg, applied in that order.

    Returns
    -------
    :
        (3, 3) matrix from the lab frame to the rotated frame.
    """
    cx, cy, cz = np.cos(np.radians([rx, ry, rz]))
    sx, sy, sz = np.sin(np.radians([rx, ry, rz]))
    return np.array(
        [
            [cy * cz, sx * sy * cz + cx * sz, sx * sz - cx * sy * cz],
            [-cy * sz, cx * cz - sx * sy * sz, sx * cz + cx * sy * sz],
            [sy, -sx * cy, cx * cy],
        ]
    )


def place_absolute(
    components: Mapping[str, "Component"], instr_param_dict: Mapping[str, float]
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Place every component in the lab frame, as McStas ABSOLUTE.

    Each component's ``at``/``rotated``/``relative`` chain is followed back to the lab frame.

    Parameters
    ----------
    components:
        Components by name, in beam order. A component placed ``relative`` to another must come
        after it.
    instr_param_dict:
        Instrument parameter values by McStas name. They resolve the expressions in ``at`` and
        ``rotated``, e.g. ``"omega"`` or ``"-chi"``.

    Returns
    -------
    :
        ``(position, rotation)`` by component name: the position in m, and the rotation matrix
        from the lab frame to the component's frame.
    """
    at_rotated_absolute_dict: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for name, c in components.items():
        parent_pos, parent_rot = (
            at_rotated_absolute_dict[c.relative] if c.relative else (np.zeros(3), np.eye(3))
        )
        at = np.array([evaluate(v, instr_param_dict) for v in c.at], float)
        rx, ry, rz = (evaluate(v, instr_param_dict) for v in c.rotated)
        rot = rotation(rx, ry, rz) @ parent_rot
        at_rotated_absolute_dict[name] = (parent_pos + parent_rot.T @ at, rot)
    return at_rotated_absolute_dict
