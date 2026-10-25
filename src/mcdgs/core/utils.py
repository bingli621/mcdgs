"""Small helpers for instruments.

- :func:`format_value` formats a value briefly, for display.
"""

from typing import Any


def format_value(value: Any) -> str:
    """Format a value briefly, for display.

    Parameters
    ----------
    value:
        A number, a list (formatted item by item), or anything else.

    Returns
    -------
    :
        A number with at most 6 significant digits, a list in brackets, or ``str(value)``.

    Examples
    --------
    >>> format_value(1 / 3)
    '0.333333'
    >>> format_value([0, 0.5, "omega"])
    '[0, 0.5, omega]'
    """
    if isinstance(value, list):
        return "[" + ", ".join(format_value(v) for v in value) + "]"
    if isinstance(value, float | int) and not isinstance(value, bool):
        return f"{value:.6g}"
    return str(value)
