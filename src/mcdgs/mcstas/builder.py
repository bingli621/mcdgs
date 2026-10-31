"""Build the McStasScript instrument of any :class:`.Instrument`.

The instrument parameters are added first, then the components one by one in beam order, since
McStas positions are relative to earlier components. Each component is:

- added as its ``component_name``, named after its attribute in the instrument, with its instance
  keywords (``at``, ``rotated``, ``relative``, ``when``, ...);
- given its McStas parameters (``Component.mcstas_parameters``: its fields plus what it
  computes), written as McStas expects them: see :func:`_formatted_value`.

A field McStas does not know is left out, since the component only computes from it. A computed
parameter McStas does not know raises.

Coordinate system: McStas, z along the beam, y up, x to the left of the beam.
"""

from dataclasses import asdict, fields
from pathlib import Path
from typing import Any

from mcstasscript.interface import instr

from mcdgs.core.base import Instrument
from mcdgs.core.components import Component
from mcdgs.core.params import InstrParam


def default_input_dir(instrument: Instrument, work_dir: str | Path | None = None) -> Path:
    """Return the folder for the McStas-generated files of an instrument.

    This is McStasScript's ``input_path``, where mcrun runs, and holds the .instr, .c and
    compiled binary. Each instrument gets its own folder so they do not overwrite each other.

    Parameters
    ----------
    instrument:
        The instrument.
    work_dir:
        Work folder. If ``None``, the current folder is used.

    Returns
    -------
    :
        ``work_dir/mcstas_input/<instrument name>``.
    """
    return Path(work_dir or Path.cwd()) / "mcstas_input" / instrument.name


def _instr_parameter(param: InstrParam) -> dict[str, Any]:
    """Return the keyword arguments of McStasScript's ``add_parameter`` for a parameter.

    McStasScript writes a number as it is, as a ``double``, but a string needs its type and C
    quotes.

    Parameters
    ----------
    param:
        The instrument parameter.

    Returns
    -------
    :
        Its ``name``, ``value``, ``unit`` and ``comment``, and ``type`` if it is a string.
    """
    if isinstance(param.value, str):
        return asdict(param) | {"type": "string", "value": f'"{param.value}"'}
    return asdict(param)


def _formatted_value(value: Any, mcstas_type: str) -> Any:
    """Write a value of a component as McStas takes it.

    Parameters
    ----------
    value:
        The value of a McStas parameter.
    mcstas_type:
        The type McStas declares for it, such as ``"string"``.

    Returns
    -------
    :
        A string in C quotes (``None`` as McStas's ``"NULL"``), a list as a C initializer in
        braces, and anything else, such as a number or an expression, unchanged.

    Examples
    --------
    >>> _formatted_value("W", "string")
    '"W"'
    >>> _formatted_value(None, "string")
    '"NULL"'
    >>> _formatted_value([0.99, 0.0219], "vector")
    '{0.99,0.0219}'
    >>> _formatted_value("-chi", "double")
    '-chi'
    """
    if mcstas_type == "string":
        return f'"{"NULL" if value is None else value}"'
    if isinstance(value, list):
        return "{" + ",".join(str(v) for v in value) + "}"
    return value


def _formatted_parameters(
    mcstas_comp: Any, component: Component, params: dict[str, Any]
) -> dict[str, Any]:
    """Return the McStas parameters of a component that its McStasScript component takes.

    Each entry of ``params`` (from ``component.mcstas_parameters``) must be one of:

    - a McStas parameter: kept, formatted by :func:`_formatted_value`;
    - a field McStas does not take: left out, as the component only computes from it (e.g. the
      crystal's ``lattice_params``, put into its ``ax`` ... ``cz``).

    Anything else was computed under a name McStas does not take, likely a typo.

    Parameters
    ----------
    mcstas_comp:
        The McStasScript component.
    component:
        The component.
    params:
        Its McStas parameters.

    Returns
    -------
    :
        The values to set, by McStas parameter name.

    Raises
    ------
    ValueError
        Naming every computed parameter McStas does not take.
    """
    mcstas_names = set(mcstas_comp.parameter_names)  # what the McStas component takes
    field_names = {f.name for f in fields(component)}
    if unknown := sorted(params.keys() - mcstas_names - field_names):
        raise ValueError(
            f"{component.name} ({component.component_name}): no McStas parameter {unknown}"
        )
    return {
        name: _formatted_value(value, mcstas_comp.parameter_types[name])
        for name, value in params.items()
        if name in mcstas_names  # else a field the component computes from: not for McStas
    }


def build_mcstas(
    instrument: Instrument, input_path: str | Path | None = None
) -> instr.McStas_instr:
    """Build the McStasScript instrument of an instrument.

    Parameters
    ----------
    instrument:
        The instrument to build.
    input_path:
        Folder for the McStas files, created if missing. If ``None``,
        :func:`default_input_dir` is used.

    Returns
    -------
    :
        The McStasScript instrument, with the instrument parameters and all components.
    """
    input_path = Path(input_path or default_input_dir(instrument))
    input_path.mkdir(parents=True, exist_ok=True)

    mcstas_instr = instr.McStas_instr(instrument.name, input_path=str(input_path))
    # a run sets its own values with `set_parameters`
    for param in instrument.instr_params:
        mcstas_instr.add_parameter(**_instr_parameter(param))

    for name, component in instrument.components.items():
        mcstas_comp = mcstas_instr.add_component(
            name,
            component.component_name,
            AT=list(component.at),
            ROTATED=list(component.rotated),
            RELATIVE=component.relative or "ABSOLUTE",  # of both AT and ROTATED
            WHEN=component.when,
            GROUP=component.group,
            JUMP=component.jump,
            SPLIT=component.split,
            EXTEND=component.extend,
        )
        params = component.mcstas_parameters()  # its fields, and what it computes
        mcstas_comp.set_parameters(_formatted_parameters(mcstas_comp, component, params))
    return mcstas_instr
