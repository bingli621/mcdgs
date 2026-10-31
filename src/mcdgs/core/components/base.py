"""Base class of every component."""

from dataclasses import dataclass, field, fields
from typing import Any, ClassVar

from mcdgs.core.utils import format_value


@dataclass(kw_only=True, repr=False)  # repr: its own, below; each subclass keeps it
class Component:
    """Base class of a McStas component.

    A subclass is named after the McStas component it stands for (``Slit``, ``DiskChopper``,
    ...). Its fields are its McStas parameters, plus the instance keywords below.

    Parameters
    ----------
    name:
        McStas instance name, and the name a ``relative`` refers to. If ``None``, it is unset.
    at:
        Position (x, y, z) in m, relative to ``relative``, as McStas AT.
    rotated:
        Rotation (about x, y, z) in deg, relative to ``relative``, as McStas ROTATED.
    relative:
        Name of the component that ``at`` and ``rotated`` are relative to.
        If ``None``, they are absolute (lab frame).
    when:
        Condition under which the component acts, a C expression (McStas WHEN, without the
        keyword or brackets). If ``None``, it always acts.
    group:
        Name of the McStas GROUP the component belongs to. If ``None``, it is in none.
    jump:
        McStas JUMP, without the keyword: a target and its condition, such as
        ``"myself ITERATE 5"``. If ``None``, there is no jump.
    split:
        Number of times each neutron is repeated on reaching the component (McStas SPLIT), or
        an expression of instrument parameters. If ``None``, it is not split.
    extend:
        C code run after the component (McStas EXTEND, without the brackets). If ``None``,
        there is none.

    Attributes
    ----------
    KEYWORDS:
        The fields that are McStas instance keywords, not McStas parameters.
    component_name:
        McStas component, the class name.

    Notes
    -----
    ``name`` is not compared, and ``component_name`` is not a field, so it is neither saved nor
    compared.

    A field typed to take a ``str`` as well as a number (``at``, ``rotated``, ``split``, and
    some of a subclass's) may be given an instrument parameter's name, or an expression of them,
    instead of a number, as in McStas. Other fields take numbers only. A component does not know
    which of its values a run sets: the instrument declares its instrument parameters.

    A McStas string parameter is a plain ``str`` (``sector="W"``), never one in C quotes: the
    McStas builder adds the quotes. A string parameter that McStas leaves unset (``"NULL"``) is
    ``None``, which the builder writes as McStas's ``"NULL"``.

    A component does not know the file it came from. A value that follows from its own values is
    computed when used; one that depends on other components is left to the instrument.

    A subclass that computes McStas parameters overrides :meth:`mcstas_parameters` and adds them
    to ``super()``'s. It computes them from its own values, never from numbers of its own. The
    fields it computes from stay in the result; the builder leaves out those McStas does not
    know.
    """

    KEYWORDS: ClassVar[tuple[str, ...]] = (
        "name",
        "at",
        "rotated",
        "relative",
        "when",
        "group",
        "jump",
        "split",
        "extend",
    )

    name: str | None = field(default=None, compare=False)
    at: list[float | str] = field(default_factory=lambda: [0, 0, 0])
    rotated: list[float | str] = field(default_factory=lambda: [0, 0, 0])
    relative: str | None = None
    when: str | None = None
    group: str | None = None
    jump: str | None = None
    split: int | str | None = None
    extend: str | None = None

    @property
    def component_name(self) -> str:
        """McStas component, the class name."""
        return type(self).__name__

    def mcstas_parameters(self) -> dict[str, Any]:
        """Return the McStas parameters of the component.

        Returns
        -------
        :
            Values by McStas parameter name: every field except the instance keywords
            (:attr:`KEYWORDS`), which McStas takes as the instance name, AT, ROTATED, ... instead.
        """
        return {f.name: getattr(self, f.name) for f in fields(self) if f.name not in self.KEYWORDS}

    def __repr__(self) -> str:
        """Show the name, the McStas component and its values as a table, one row per field.

        The placement is always shown; the other instance keywords only when they are set.
        """
        from tabulate import tabulate  # only when printed: importing does not need it

        optional = set(self.KEYWORDS) - {"at", "rotated", "relative"}
        title = self.component_name if self.name is None else f"{self.name} ({self.component_name})"
        rows = [
            (f.name, format_value(getattr(self, f.name)))
            for f in fields(self)
            if not (f.name in optional and getattr(self, f.name) is None) and f.name != "name"
        ]
        return f"{title}:\n{tabulate(rows, ('parameter', 'value'), disable_numparse=True)}"
