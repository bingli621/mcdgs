"""Build the McStasScript instrument of any :class:`.Instrument`.

The instrument parameters are added first, then the components one by one in beam order, since
McStas positions are relative to earlier components. Each component is:

- added as its ``component_name``, named after its attribute in the instrument, with its instance
  keywords (``at``, ``rotated``, ``relative``, ``when``, ...);
- given its McStas parameters (``Component.mcstas_parameters``: its fields plus what it
  computes), written as McStas expects them: see :func:`_parameter_value`. A list, the value
  of a McStas ``vector``, is declared as an array named ``<component>_<parameter>`` in DECLARE
  and passed by that name. McStas 3.10.2 sizes a vector written in place by its first instance,
  so the edges of a second chopper with more slits would overflow.

The chopper train of a ``Polygon_ESS_butterfly`` is an array of chopper-lib's
``chopper_parameters``, which McStas takes as a ``vector``: see :func:`_declare_chopper_train`.

A field McStas does not know is left out, since the component only computes from it. A computed
parameter McStas does not know raises.

Coordinate system: McStas, z along the beam, y up, x to the left of the beam.
"""

from dataclasses import asdict, fields
from pathlib import Path
from typing import Any

from mcstasscript.interface import instr

from mcdgs.core import InstrParam, Instrument
from mcdgs.core.components import Component, Polygon_ESS_butterfly
from mcdgs.mcstas.workspace import Workspace


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


def _declare_chopper_train(
    mcstas_instr: instr.McStas_instr, component: Polygon_ESS_butterfly
) -> str | int:
    """Declare the chopper train of a ``Polygon_ESS_butterfly``, and return what McStas takes.

    McStas takes the train as a ``vector``, a ``double *``: the address of an array of
    chopper-lib's ``chopper_parameters``. DECLARE holds the edges of each disk, the array
    ``<component>_choppers`` and the pointer ``<component>_choppers_pointer``. INITIALIZE fills
    the array, since a disk's ``speed`` and ``delay`` may be instrument parameters, and points
    the pointer at it: McStas runs it before the INITIALIZE of the components.

    Both are C code inside ``#ifdef CHOPPER_LIB_VERSION``, which the SHARE of the source
    defines, before DECLARE: a run from ``sample_pos`` leaves the source out, and with it the
    type ``chopper_parameters``, so the code drops out with it.

    Parameters
    ----------
    mcstas_instr:
        The McStasScript instrument, to declare the arrays in.
    component:
        The source.

    Returns
    -------
    :
        The name of the pointer, or without a disk 0, the null pointer: McStasScript refuses
        ``NULL``, a name it does not know.
    """
    if not component.choppers:
        return 0
    train = f"{component.name}_choppers"
    pointer = f"{train}_pointer"
    declare = [f"chopper_parameters {train}[{len(component.choppers)}];"]
    initialize = []
    for i, disk in enumerate(component.choppers):
        edges = f"{train}_edges_{i}"
        values = ", ".join(str(edge) for edge in disk["edges"])
        declare.append(f"double {edges}[{len(disk['edges'])}] = {{{values}}};")
        initialize.append(
            f"{train}[{i}] = (chopper_parameters){{.speed = {disk['speed']}, "
            f".delay = {disk['delay']}, .beam = {disk['beam']}, "
            f".edge_count = {len(disk['edges'])}, .edges = {edges}, "
            f".path = {disk['path']}, .aperture = {disk['aperture']}}};"
        )
    declare.append(f"double * {pointer};")
    initialize.append(f"{pointer} = (double *) {train};")
    mcstas_instr.append_declare("\n".join(["#ifdef CHOPPER_LIB_VERSION", *declare, "#endif"]))
    mcstas_instr.append_initialize("\n".join(["#ifdef CHOPPER_LIB_VERSION", *initialize, "#endif"]))
    return pointer


def _parameter_value(
    mcstas_instr: instr.McStas_instr,
    mcstas_comp: Any,
    component: Component,
    name: str,
    value: Any,
) -> Any:
    """Return the value of one McStas parameter as its McStasScript component takes it.

    Parameters
    ----------
    mcstas_instr:
        The McStasScript instrument, to declare the arrays in.
    mcstas_comp:
        The McStasScript component.
    component:
        The component.
    name:
        The McStas parameter.
    value:
        Its value.

    Returns
    -------
    :
        One of:

        - a list: the array ``<component>_<parameter>`` it is declared as, passed by that name
          (for the chopper train of a ``Polygon_ESS_butterfly``, what
          :func:`_declare_chopper_train` returns);
        - a value of a McStas ``string``: in C quotes, with ``None`` as McStas's ``"NULL"``;
        - anything else, such as a number or an expression: unchanged.
    """
    if isinstance(component, Polygon_ESS_butterfly) and name == "choppers":
        return _declare_chopper_train(mcstas_instr, component)
    if isinstance(value, list):  # a vector: declared, not written in place
        return mcstas_instr.add_declare_var(
            "double", f"{component.name}_{name}", array=len(value), value=value
        )
    if mcstas_comp.parameter_types[name] == "string":
        return f'"{"NULL" if value is None else value}"'
    return value


def _formatted_parameters(
    mcstas_instr: instr.McStas_instr,
    mcstas_comp: Any,
    component: Component,
    params: dict[str, Any],
) -> dict[str, Any]:
    """Return the McStas parameters of a component that its McStasScript component takes.

    Each entry of ``params`` (from ``component.mcstas_parameters``) must be one of:

    - a McStas parameter: kept, as :func:`_parameter_value` writes it;
    - a field McStas does not take: left out, as the component only computes from it (e.g. the
      crystal's ``lattice_params``, put into its ``ax`` ... ``cz``).

    Anything else was computed under a name McStas does not take, likely a typo.

    Parameters
    ----------
    mcstas_instr:
        The McStasScript instrument, to declare the arrays in.
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
        name: _parameter_value(mcstas_instr, mcstas_comp, component, name, value)
        for name, value in params.items()
        if name in mcstas_names  # not a field the component only computes from
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
        Folder for the McStas files, created if missing. If ``None``, the build folder of the
        instrument in :meth:`.Workspace.find`.

    Returns
    -------
    :
        The McStasScript instrument, with the instrument parameters and all components.
    """
    input_path = Path(input_path or Workspace.find().build_folder(instrument.name))
    input_path.mkdir(parents=True, exist_ok=True)

    mcstas_instr = instr.McStas_instr(instrument.name, input_path=str(input_path))
    # a run sets its own values with `set_parameters`
    for param in instrument.instr_params:
        mcstas_instr.add_parameter(**_instr_parameter(param))

    for name, component in instrument.components.items():  # in beam order
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
        mcstas_comp.set_parameters(
            _formatted_parameters(mcstas_instr, mcstas_comp, component, params)
        )
    return mcstas_instr
