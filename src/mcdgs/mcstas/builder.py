"""Build the McStasScript instrument of any :class:`.Instrument`.

The instrument parameters are added first, then the components one by one in beam order, since
McStas positions are relative to earlier components. Each component is:

- added as its ``component_name``, named after its attribute in the instrument, with its instance
  keywords (``at``, ``rotated``, ``relative``, ``when``, ...);
- given its McStas parameters (``Component.mcstas_parameters``: its fields plus what it
  computes), written as McStas expects them: see :func:`_formatted_value`. A list, the value
  of a McStas ``vector``, is declared as an array named ``<component>_<parameter>`` in DECLARE
  and passed by that name. McStas 3.10.2 sizes a vector written in place by its first instance,
  so the edges of a second chopper with more slits would overflow.

A field McStas does not know is left out, since the component only computes from it. A computed
parameter McStas does not know raises.

Coordinate system: McStas, z along the beam, y up, x to the left of the beam.
"""

import warnings
from dataclasses import asdict, fields
from pathlib import Path
from typing import Any

from mcstasscript.helper.component_reader import ComponentReader
from mcstasscript.interface import instr

from mcdgs.core.base import Instrument
from mcdgs.core.components import Component, ESS_butterfly, Monitor_nD, NXdisk_chopper
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
        A string in C quotes (``None`` as McStas's ``"NULL"``), and anything else, such as a
        number or an expression, unchanged. A list is declared instead: see
        :func:`_formatted_parameters`.

    Examples
    --------
    >>> _formatted_value("W", "string")
    '"W"'
    >>> _formatted_value(None, "string")
    '"NULL"'
    >>> _formatted_value("-chi", "double")
    '-chi'
    """
    if mcstas_type == "string":
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

    - a McStas parameter: kept, formatted by :func:`_formatted_value`, or, if it is a list,
      declared as the array ``<component>_<parameter>`` and passed by that name;
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
    values = {}
    for name, value in params.items():
        if name not in mcstas_names:  # a field the component computes from: not for McStas
            continue
        if isinstance(value, list):  # a vector: declared, not written in place
            values[name] = mcstas_instr.add_declare_var(
                "double", f"{component.name}_{name}", array=len(value), value=value
            )
        else:
            values[name] = _formatted_value(value, mcstas_comp.parameter_types[name])
    return values


def build_mcstas(
    instrument: Instrument, input_path: str | Path | None = None, *, tof_train: bool = False
) -> instr.McStas_instr:
    """Build the McStasScript instrument of an instrument.

    Parameters
    ----------
    instrument:
        The instrument to build.
    input_path:
        Folder for the McStas files, created if missing. If ``None``,
        :func:`default_input_dir` is used.
    tof_train:
        Opt into the experimental T-Rex TOF train component overrides. Generated files and
        overrides are loaded from the bundled component directory; generated files remain in
        ``input_path``. OpenACC and OpenMP are unsupported. MPI ranks use
        independent adaptive feedback. Train samples are correlated; monitor error estimates
        must not be interpreted as independent-ray uncertainties.

    Returns
    -------
    :
        The McStasScript instrument, with the instrument parameters and all components.
    """
    input_path = Path(input_path or default_input_dir(instrument))
    input_path.mkdir(parents=True, exist_ok=True)

    if tof_train:
        from mcdgs.instruments.trex.model import TRex

        if type(instrument) is not TRex:
            raise ValueError("TOF trains are restricted to the reviewed T-Rex model")
        reference = TRex(params=instrument.params).components
        actual = instrument.components
        if list(actual) != list(reference) or any(
            type(actual[name]) is not type(reference[name]) for name in reference if name in actual
        ):
            raise ValueError("TOF trains require the reviewed sample-free T-Rex component sequence")
        if any(c.extend for c in actual.values()):
            raise ValueError("TOF trains do not support custom EXTEND code")
        if any(actual[name].when != reference[name].when for name in reference):
            raise ValueError("TOF trains do not support custom WHEN conditions")
        monitor = actual["pre_sample_monitor"]
        expected_monitor = reference["pre_sample_monitor"]
        if (monitor.at, monitor.rotated, monitor.relative, monitor.when) != (
            expected_monitor.at,
            expected_monitor.rotated,
            expected_monitor.relative,
            expected_monitor.when,
        ):
            raise ValueError("TOF trains require the reviewed pre-sample monitor placement")
        if any(
            getattr(monitor, name) != getattr(expected_monitor, name)
            for name in ("xwidth", "yheight", "restore_neutron")
        ):
            raise ValueError("TOF trains require the reviewed pre-sample monitor geometry")
        components_dir = (
            Path(__file__).resolve().parents[1] / "instruments/trex/tof_train_components"
        )
        required = (
            "ESS_butterfly.comp",
            "ESS_butterfly-lib.c",
            "ESS_butterfly-lib.h",
            "NXdisk_chopper.comp",
            "chopper-lib.c",
            "chopper-lib.h",
            "Monitor_nD.comp",
            "trex-tof-train-monitor_nd-lib.h",
            "trex-tof-train-monitor_nd-lib.c",
            "trex-tof-train-geometry-mesh-lib.h",
            "trex-tof-train-geometry-mesh-lib.c",
            "trex-tof-train-butterfly-geometry.c",
        )
        missing = [name for name in required if not (components_dir / name).is_file()]
        if missing:
            raise FileNotFoundError(f"TOF train component overrides missing: {missing}")
        collisions = [name for name in required if (input_path / name).exists()]
        if collisions:
            raise FileExistsError(
                f"TOF train search may be shadowed by files in {input_path}: {collisions}"
            )
    mcstas_instr = instr.McStas_instr(instrument.name, input_path=str(input_path))
    if tof_train:
        # McStasScript's instrument input_path also controls generated files and its beam
        # database, so use the override directory only for component discovery.
        mcstas_instr.component_reader = ComponentReader(
            mcstas_instr._run_settings["package_path"], input_path=str(components_dir)
        )
        for name in ("ESS_butterfly", "NXdisk_chopper", "Monitor_nD"):
            resolved = Path(mcstas_instr.component_reader.component_path[name]).resolve()
            if resolved != (components_dir / f"{name}.comp").resolve():
                raise RuntimeError(f"TOF train {name} override was not loaded: {resolved}")
        mcstas_instr.add_search(str(components_dir))
        components = list(instrument.components.items())
        sources = [i for i, (_, c) in enumerate(components) if isinstance(c, ESS_butterfly)]
        choppers = [i for i, (_, c) in enumerate(components) if isinstance(c, NXdisk_chopper)]
        if len(sources) != 1 or not choppers or sources[0] >= choppers[0]:
            raise ValueError("TOF trains require one ESS source followed by disk choppers")
        source = components[sources[0]][1]
        if source.when:
            raise ValueError("TOF trains require an unconditional source")
        if isinstance(source, ESS_butterfly) and source.tfocus_width != 0:
            raise ValueError("TOF trains do not support focused emission time")
        targets = [
            name
            for name, c in components[choppers[-1] + 1 :]
            if isinstance(c, Monitor_nD) and not c.when
        ]
        if not targets:
            raise ValueError("TOF trains require a monitor after the last chopper")
        target = "pre_sample_monitor"
        if target not in targets:
            raise ValueError("TOF trains require the pre-sample monitor after the last chopper")
        mcstas_instr.add_parameter("N_trains_par", value=300, type="int")
        mcstas_instr.add_parameter("target_tsplit", value=10)
        for var_type, name in (
            ("double *", "t_offset"),
            ("double *", "p_trains"),
            ("double", "p_last_time_manipulation"),
            ("int", "adaptive_N"),
            ("int", "N_active"),
            ("long *", "train_total_arrived"),
            ("long *", "train_total_N_sent"),
            ("long *", "train_total_rays_sent"),
        ):
            mcstas_instr.add_user_var(var_type, name)
        if any(c.jump for _, c in components):
            raise ValueError("TOF trains do not support jumps")
        if any(c.split for _, c in components):
            warnings.warn(
                "TOF train SPLIT must be after the last chopper, with no absorbing "
                "components that modify the train in downstream branches; shared train "
                "arrays and time-dependent components may otherwise give incorrect results.",
                UserWarning,
                stacklevel=2,
            )
        mcstas_instr.add_declare_var("double *", "tof_train_offsets")
        mcstas_instr.add_declare_var("double *", "tof_train_weights")
        for name in ("tof_train_arrived", "tof_train_sent", "tof_train_rays"):
            mcstas_instr.add_declare_var("long", name, value=0)
        mcstas_instr.add_declare_var("int", "N_trains")
        mcstas_instr.append_initialize("""#if defined(OPENACC) || defined(_OPENACC) || defined(_OPENMP)
#error TOF train integration requires serial particle tracing
#endif""")
        mcstas_instr.append_initialize("""N_trains = N_trains_par;
if (N_trains < 1) { fprintf(stderr, "N_trains must be positive\\n"); exit(EXIT_FAILURE); }
tof_train_offsets = malloc(N_trains*sizeof(double));
tof_train_weights = malloc(N_trains*sizeof(double));
if (!tof_train_offsets || !tof_train_weights) exit(EXIT_FAILURE);""")
        mcstas_instr.append_finally("free(tof_train_offsets); free(tof_train_weights);")
    # a run sets its own values with `set_parameters`
    for param in instrument.instr_params:
        mcstas_instr.add_parameter(**_instr_parameter(param))

    for name, component in instrument.components.items():
        if tof_train and isinstance(component, ESS_butterfly):
            mcstas_instr.add_component(
                "tof_train_init",
                "Arm",
                AT=[0, 0, 0],
                RELATIVE="ABSOLUTE",
                EXTEND="""t_offset = tof_train_offsets;
p_trains = tof_train_weights;
p_last_time_manipulation = 0;
adaptive_N = N_trains;
N_active = 0;
train_total_arrived = &tof_train_arrived;
train_total_N_sent = &tof_train_sent;
train_total_rays_sent = &tof_train_rays;""",
            )
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
        if tof_train and isinstance(component, ESS_butterfly):
            if "target_tsplit" not in mcstas_comp.parameter_names:
                raise RuntimeError("TOF train source override was not loaded")
            params["target_tsplit"] = "target_tsplit"
        if tof_train and isinstance(component, Monitor_nD):
            if not {"tsplit", "adaptive_target"} <= set(mcstas_comp.parameter_names):
                raise RuntimeError("TOF train monitor override was not loaded")
            params["adaptive_target"] = int(name == target)
        mcstas_comp.set_parameters(
            _formatted_parameters(mcstas_instr, mcstas_comp, component, params)
        )
    return mcstas_instr
