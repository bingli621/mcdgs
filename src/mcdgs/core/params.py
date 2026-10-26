"""Parameter types shared by all instruments.

- :func:`param` declares a field of a control-parameters dataclass, with a unit and a comment
  that io/ writes next to its value. Every instrument's ``params_class`` is built from them.
- :class:`InstrParam` is one McStas instrument parameter, which a run may set. An instrument's
  ``instr_params`` is a list of them, computed from its control parameters.

Components stay general: a component value that an instrument parameter provides is given as
the parameter's name, or an expression of them, as in McStas (``rotated: [0, "omega", 0]``).

Design parameters are the components' own values. How they are stored in a file belongs to
io/, and how McStas runs belongs to mcstas/.
"""

from dataclasses import dataclass, field
from typing import Any


@dataclass(frozen=True)
class InstrParam:
    """A McStas instrument parameter.

    Parameters
    ----------
    name:
        McStas name.
    value:
        Default value.
    unit:
        Unit, for display.
    comment:
        Description, for display.
    """

    name: str
    value: object = None
    unit: str = ""
    comment: str = ""


def param(default: Any, unit: str, comment: str) -> Any:
    """Declare a field of a ``params_class``.

    Parameters
    ----------
    default:
        Default value.
    unit:
        Unit, stored in the field's metadata.
    comment:
        Description, stored in the field's metadata.

    Returns
    -------
    :
        A dataclass field to assign in the class body. It is typed ``Any``, as
        ``dataclasses.field`` is, so that a type checker accepts it as the field's default.

    Examples
    --------
    >>> from dataclasses import dataclass, fields
    >>> @dataclass(frozen=True)
    ... class Params:
    ...     ei: float = param(25.0, "meV", "incident energy")
    >>> fields(Params)[0].metadata["unit"]
    'meV'
    """
    return field(default=default, metadata={"unit": unit, "comment": comment})
