"""The ISCS frame of the CSPEC specifications, in which every guide and gap is placed.

The specifications (CSPEC nl-specs 1.E, the .xlsx file in ``reference/``) give each element of
the guide in ISCS, the frame of the insert of beamport W3, in mm: x along the insert, y across
it, z up. ``arm_iscs`` is that frame in the McStas frame of the source, whose origin is the
focal point of the sector and whose z is along beamline W3. The insert is 0.5 deg off the
beamline, toward -x.

Each element is given by its start (x, y, z) and two angles, chi and phi, in deg: the trajectory
of the guide turns down and back up by chi (the vertical S-bender) and sideways by phi (the
horizontal curve). The angles are those of the trajectory at the middle of the element, and the
sheet takes the element straight from its start at them. :func:`guide` and :func:`arm` place
it so, relative to ``arm_iscs``.

In ``arm_iscs``, McStas x is ISCS y, McStas y is ISCS z and McStas z is ISCS x, which keeps the
frame right-handed. The sheet's element runs along (cos phi cos chi, -sin phi, sin chi) in ISCS,
so McStas ROTATED (-chi, -phi, 0) is its direction, to a factor cos phi on the vertical.
"""

from typing import Any

from mcdgs.core.components import Arm, Guide_gravity

BEAMLINE_ANGLE = 138.0  # deg, cell K3 of Coords-Trans, beamline W3
ISCS_ANGLE = 137.5  # deg, cell T3 of Coords-Trans, the insert of beamport W3
ISCS_ORIGIN_X = 0.03271038285817347  # m, cell N10 of Coords-Trans
ISCS_ORIGIN_Z = 0.02437274816822491  # m, cell P10 of Coords-Trans

MM_PER_M = 1000  # the sheet's lengths are in mm


def arm_iscs() -> Arm:
    """Create the arm of the ISCS frame, which every guide and gap hangs off.

    Returns
    -------
    :
        A new arm, at the ISCS origin, rotated from the beamline to the insert.
    """
    return Arm(
        at=[ISCS_ORIGIN_X, 0, ISCS_ORIGIN_Z],
        rotated=[0, ISCS_ANGLE - BEAMLINE_ANGLE, 0],
    )


def _placement(x: float, y: float, z: float, chi: float, phi: float) -> dict[str, Any]:
    """Return the McStas placement of a point and direction of the sheet, in ``arm_iscs``.

    Parameters
    ----------
    x, y, z:
        Position in ISCS, in mm.
    chi, phi:
        Vertical and horizontal angle of the trajectory, in deg.

    Returns
    -------
    :
        ``at``, ``rotated`` and ``relative`` of the component.
    """
    return {
        "at": [y / MM_PER_M, z / MM_PER_M, x / MM_PER_M],
        "rotated": [-chi, -phi, 0],
        "relative": "arm_iscs",
    }


def arm(x: float, y: float, z: float, chi: float = 0.0, phi: float = 0.0) -> Arm:
    """Create an arm at a point of the sheet, along the trajectory there.

    Parameters
    ----------
    x, y, z:
        Position in ISCS, in mm.
    chi, phi:
        Vertical and horizontal angle of the trajectory, in deg.

    Returns
    -------
    :
        A new arm, relative to ``arm_iscs``.
    """
    return Arm(**_placement(x, y, z, chi, phi))


def guide(
    *,
    x: float,
    y: float,
    z: float,
    chi: float,
    phi: float,
    length: float,
    width_in: float,
    width_out: float,
    height_in: float,
    height_out: float,
    mtop: float,
    mbottom: float,
    mleft: float,
    mright: float,
) -> Guide_gravity:
    """Create a guide element of the sheet, in mm and deg as the sheet gives it.

    Parameters
    ----------
    x, y, z:
        Start in ISCS, in mm.
    chi, phi:
        Vertical and horizontal angle of the trajectory at the middle of the element, in deg.
    length:
        Length in mm.
    width_in, width_out, height_in, height_out:
        Width and height at the entry and the exit, in mm.
    mtop, mbottom, mleft, mright:
        m-value of each mirror.

    Returns
    -------
    :
        A new guide, relative to ``arm_iscs``.
    """
    return Guide_gravity(
        w1=width_in / MM_PER_M,
        h1=height_in / MM_PER_M,
        w2=width_out / MM_PER_M,
        h2=height_out / MM_PER_M,
        l=length / MM_PER_M,
        mtop=mtop,
        mbottom=mbottom,
        mleft=mleft,
        mright=mright,
        **_placement(x, y, z, chi, phi),
    )
