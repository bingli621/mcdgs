"""Small helpers for instruments.

- :func:`format_value` formats a value briefly, for display.
- :func:`evaluate` computes the number a component value stands for. A McStas expression of
  instrument parameters (``"-chi"``, ``"2.0*omega + 5"``) is parsed with ``ast`` and computed as
  McStas would. Only what McStas offers is allowed, and no code is run.
"""

import ast
import math
import operator
from collections.abc import Callable, Mapping
from typing import Any

# what an expression may use: McStas's constants and its C math functions
CONSTANTS: dict[str, float] = {"PI": math.pi, "DEG2RAD": math.pi / 180, "RAD2DEG": 180 / math.pi}
FUNCTIONS: dict[str, Callable[..., float]] = {
    "sin": math.sin,
    "cos": math.cos,
    "tan": math.tan,
    "asin": math.asin,
    "acos": math.acos,
    "atan": math.atan,
    "atan2": math.atan2,
    "sqrt": math.sqrt,
    "fabs": math.fabs,
    "exp": math.exp,
    "log": math.log,
    "pow": math.pow,
}
_BINARY: dict[type[ast.operator], Callable[[float, float], float]] = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Mult: operator.mul,
    ast.Div: operator.truediv,
}
_UNARY: dict[type[ast.unaryop], Callable[[float], float]] = {
    ast.USub: operator.neg,
    ast.UAdd: operator.pos,
}


def evaluate(entry: float | str, params: Mapping[str, float]) -> float:
    """Compute the number a component value stands for.

    Parameters
    ----------
    entry:
        A number, returned unchanged, or a McStas expression. An expression may use numbers,
        instrument parameters, ``CONSTANTS``, ``+ - * /``, signs, brackets and ``FUNCTIONS``.
    params:
        Instrument parameter values by name.

    Returns
    -------
    :
        The value of ``entry``.

    Raises
    ------
    ValueError
        If ``entry`` is not a valid expression, uses anything not listed above, or divides two
        integers. C would round such a division down, so write ``1.0/2``.

    Examples
    --------
    >>> params = {"omega": 30, "chi": 5}
    >>> evaluate(1.5, params)
    1.5
    >>> evaluate("omega", params)
    30
    >>> evaluate("-chi", params)
    -5
    >>> evaluate("2.0*omega + chi", params)
    65.0
    >>> evaluate("atan2(1, 1)*RAD2DEG", params)
    45.0
    """
    if not isinstance(entry, str):
        return entry

    def value(node: ast.expr) -> float:
        match node:
            case ast.Constant(value=int() | float() as number) if not isinstance(number, bool):
                return number
            case ast.Name(id=name) if name in params:
                return params[name]
            case ast.Name(id=name) if name in CONSTANTS:
                return CONSTANTS[name]
            case ast.BinOp(left=left, op=op, right=right) if type(op) in _BINARY:
                a, b = value(left), value(right)
                if isinstance(op, ast.Div) and isinstance(a, int) and isinstance(b, int):
                    msg = f"{entry!r}: C divides integers rounding down; write 1.0/2"
                    raise ValueError(msg)  # noqa: TRY004 (a bad expression, like the others)
                return _BINARY[type(op)](a, b)
            case ast.UnaryOp(op=op, operand=operand) if type(op) in _UNARY:
                return _UNARY[type(op)](value(operand))
            case ast.Call(func=ast.Name(id=name), args=args, keywords=[]) if name in FUNCTIONS:
                return FUNCTIONS[name](*(value(a) for a in args))
        raise ValueError(
            f"{entry!r}: {ast.unparse(node)!r} is not a number, an instrument parameter, a "
            "constant, + - * /, or a McStas math function"
        )

    try:
        tree = ast.parse(entry, mode="eval")
    except SyntaxError as error:
        raise ValueError(f"{entry!r} is not an expression: {error.msg}") from None
    return value(tree.body)


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
