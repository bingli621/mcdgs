"""Everything that talks to McStas: building the instrument and running simulations.

The names below are imported from here, ``from mcdgs.mcstas import Simulation``.
"""

from .builder import build_mcstas
from .mcstasscript_diagnostics import McStasScriptDiagnostics
from .nexus_io import merge_entries, numeric_parameters, read_events, read_pixels
from .reader import read_mcstas_map
from .runs import Run, RunType, monitor_file, recorded_instr_params
from .scan import compile_instrument, random_seed, run_points_in_parallel
from .simulation import Simulation
from .workspace import Workspace

__all__ = [
    "McStasScriptDiagnostics",
    "Run",
    "RunType",
    "Simulation",
    "Workspace",
    "build_mcstas",
    "compile_instrument",
    "merge_entries",
    "monitor_file",
    "numeric_parameters",
    "random_seed",
    "read_events",
    "read_mcstas_map",
    "read_pixels",
    "recorded_instr_params",
    "run_points_in_parallel",
]
