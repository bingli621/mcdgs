"""Simulate an instrument into McStas NeXus (HDF5) output, read with :func:`.read_pixels`.

A :class:`Simulation` holds the instrument to simulate and the :class:`.Workspace` its runs
write into. It offers:

- ``run``: one run, source -> detector;
- ``scan``: one run per value of an instrument parameter, source -> detector, in parallel;
- ``run_to_sample``: source -> ``sample_pos`` once, the beam dumped into an MCPL beam file;
- ``run_from_sample``: one run from that beam file, ``sample_pos`` -> detector;
- ``scan_from_sample``: one run per value from that beam file, in parallel.

Everything runs through McStasScript's public API (``settings``, ``set_parameters``,
``run_to``, ``run_from``, ``backengine``). A scan runs ``backengine`` once per value and merges
the outputs into one file, after compiling the instrument once: see :mod:`.mcstas.scan`.

Each of them returns the :class:`.Run` it made, which is also left as ``run.json`` in the
folder of the run.

``prepare`` is the step they all share, for a run that McStasScript starts itself, such as
``show_diagram(analysis=True)`` or the diagnostics of :class:`.McStasScriptDiagnostics`: it gives
the run its folder, its settings and its record, which the caller finishes once McStas has run.
"""

from collections.abc import Mapping
from datetime import datetime
from typing import Any

import numpy as np
from mcstasscript.interface.instr import McStas_instr
from numpy.typing import ArrayLike

from mcdgs.core import Instrument
from mcdgs.mcstas.builder import build_mcstas
from mcdgs.mcstas.runs import BEAM_FILE, Run, RunType, recorded_instr_params
from mcdgs.mcstas.scan import random_seed, run_points_in_parallel
from mcdgs.mcstas.workspace import Workspace

SAMPLE_POS = "sample_pos"  # the component a split run cuts the beam at


class Simulation:
    """Run McStas simulations of an instrument.

    Every run builds the instrument anew, into its build folder, so a change to the instrument
    applies to the next run.

    Parameters
    ----------
    instrument:
        The instrument to simulate.
    workspace:
        The folders to write into: the McStas files into its build folder of the instrument,
        each run into a folder of its runs folder. If ``None``, :meth:`.Workspace.find`, fixed
        when the simulation is created.

    Attributes
    ----------
    beam_run:
        The name of the last :meth:`run_to_sample` of this simulation, whose beam the runs from
        the sample start from unless they name another. ``None`` before the first one: a beam
        of an earlier session is never guessed, it is named.

    Notes
    -----
    Each run writes into a new folder of the runs folder, named by its start time and type:
    ``<time>_run``, ``<time>_scan_<variable>``, ``<time>_run_to_sample``,
    ``<time>_run_from_sample`` and ``<time>_scan_<variable>_from_sample``, with ``<time>`` e.g.
    ``20261022-143005``. A ``folder_name`` names the folder instead, and it is then replaced on
    every run. :meth:`.Workspace.list_runs` and :meth:`.Workspace.latest_run` find finished
    runs, and :meth:`.Workspace.clean_runs` removes older ones.

    The ``instr_params`` of a run are its instrument parameters by McStas name
    (``Instrument.instr_params``). Those not given keep the instrument's values.

    A split run (:meth:`run_to_sample`, ...) needs a component named ``sample_pos``. Its second
    part puts ``sample_pos`` at the origin, unrotated, and the beam file holds the neutrons in the
    frame of ``sample_pos``. So the frame of that part is the frame of ``sample_pos``, and a
    component after it must be placed ``relative`` to it, not in the lab frame.
    """

    def __init__(self, instrument: Instrument, workspace: Workspace | None = None) -> None:
        self.instrument = instrument
        self.workspace = workspace or Workspace.find()
        self.beam_run: str | None = None

    def prepare(
        self,
        folder_name: str | None,
        run_type: RunType = RunType.RUN,
        instr_params: Mapping[str, float] | None = None,
        *,
        ncount: float | None = None,
        seed: int | None = None,
        variable: str | None = None,
        values: list[float] | None = None,
        beam_run: str | None = None,
        nexus: bool = True,
        **settings: object,
    ) -> tuple[Run, McStas_instr]:
        """Pick the run folder and build the instrument with the settings every run shares.

        Every method of the simulation starts with it. Call it for a run that McStasScript
        starts itself, such as ``inst.show_diagram(analysis=True)``, which runs ``inst`` with
        its settings, then save the record with ``run.finish()``.

        Parameters
        ----------
        folder_name:
            Run folder name. If ``None``, the start time and what the run does, e.g.
            ``<time>_scan_omega``.
        run_type:
            What the run does.
        instr_params:
            The instrument parameters given to the run, of a scan those other than
            ``variable``. They are recorded, and set unless the run is a scan, whose points
            each set their own. If ``None``, none: all keep the instrument's values.
        ncount:
            Number of neutrons. If ``None``, as from the sample, it is not set.
        seed:
            Random seed. If ``None``, McStas picks one.
        variable, values:
            The instrument parameter of a scan and its values.
        beam_run:
            The name of the run to the sample a run from the sample starts from.
        nexus:
            True writes McStas NeXus, the ``mccode.h5`` of the record. False writes McStas's own
            format, which McStasScript's diagnostics and diagram read, and ``ms.load_data``:
            the record then names a ``data_file`` that does not exist.
        **settings:
            The run's other McStasScript settings, e.g. ``mpi``.

        Returns
        -------
        :
            The record of the run, not saved yet, and the McStasScript instrument writing into
            its folder. A scan's points write into folders of their own instead
            (:func:`.run_points_in_parallel`).
        """
        instr_params = {} if instr_params is None else instr_params
        started = datetime.now().astimezone()
        label = (
            run_type.value
            if variable is None
            else run_type.value.replace("scan", f"scan_{variable}")
        )
        folder = self.workspace.new_run_folder(folder_name, label)
        input_path = self.workspace.build_folder(self.instrument.name)
        inst = build_mcstas(self.instrument, input_path=input_path)
        if ncount is not None:
            settings["ncount"] = ncount
        inst.settings(
            NeXus=nexus,
            increment_folder_name=False,
            suppress_output=True,
            output_path=str(folder),
            force_compile=True,
            seed=seed,
            **settings,  # pyright: ignore[reportArgumentType]
        )
        if variable is None:
            inst.set_parameters(**instr_params)  # the rest keep the instrument's values
        run = Run(
            folder=folder,
            instrument=self.instrument.name,
            run_type=run_type,
            started=started,
            instr_params=recorded_instr_params(
                self.instrument.instr_params, instr_params, variable
            ),
            ncount=ncount,
            seed=seed,
            variable=variable,
            values=values,
            beam_run=beam_run,
        )
        return run, inst

    def _prepare_from_sample(
        self,
        folder_name: str | None,
        run_type: RunType,
        instr_params: Mapping[str, float],
        beam_run: str | None,
        **prepare: Any,
    ) -> tuple[Run, McStas_instr]:
        """Prepare a run from the sample: as :meth:`prepare`, cut to start from a beam file.

        Parameters
        ----------
        folder_name, run_type, instr_params:
            As in :meth:`prepare`.
        beam_run:
            Name of the run to the sample whose beam to start from. If ``None``,
            :attr:`beam_run`.
        **prepare:
            The other keywords of :meth:`prepare`.

        Returns
        -------
        :
            The record of the run, not saved yet, and the McStasScript instrument starting at
            ``sample_pos`` from the beam file.

        Raises
        ------
        ValueError
            If the instrument has no ``sample_pos``, before anything runs.
        FileNotFoundError
            If there is no run to the sample to start from, before anything runs.
        """
        self._check_sample_pos()  # before anything
        beam = self.find_beam_run(beam_run, folder_name)
        run, inst = self.prepare(folder_name, run_type, instr_params, beam_run=beam.name, **prepare)
        # in C quotes, and absolute, because mcrun runs in the build folder
        inst.run_from(SAMPLE_POS, filename=f'"{(beam.folder / BEAM_FILE).resolve()}"')
        return run, inst

    @staticmethod
    def _run_once(run: Run, inst: McStas_instr) -> Run:
        """Run McStas once, then save and return the record of the finished run."""
        inst.backengine()  # compiles in its build folder
        return run.finish()

    @staticmethod
    def _run_points(
        run: Run, inst: McStas_instr, points: list[dict[str, float]], n_jobs: int | None
    ) -> Run:
        """Run McStas once per point of a scan, then save and return the finished run's record."""
        run_points_in_parallel(inst, run.folder, points, seed=run.seed, n_jobs=n_jobs)
        return run.finish()

    def _check_sample_pos(self) -> None:
        """Check that the instrument has the component a split run cuts the beam at.

        Raises
        ------
        ValueError
            If the instrument has no ``sample_pos`` component.
        """
        if SAMPLE_POS not in self.instrument.components:
            raise ValueError(f"{self.instrument.name} has no {SAMPLE_POS} component to split at")

    def find_beam_run(self, beam_run: str | None = None, folder_name: str | None = None) -> Run:
        """Return the run to the sample that a run from the sample starts from.

        Parameters
        ----------
        beam_run:
            The name of the run to the sample, in the runs folder of the workspace. If
            ``None``, :attr:`beam_run`.
        folder_name:
            The folder name of the run from the sample, which must not be that of the run to
            the sample: the run would remove the beam it starts from.

        Returns
        -------
        :
            The run to the sample.

        Raises
        ------
        FileNotFoundError
            If both ``beam_run`` and :attr:`beam_run` are ``None``, or the run is not a
            finished run of the workspace.
        ValueError
            If the run is not a ``run_to_sample``, or ``folder_name`` is its name.
        """
        beam_run = beam_run or self.beam_run
        if beam_run is None:
            raise FileNotFoundError(
                "no run to the sample: run run_to_sample first, or name an earlier one as beam_run"
            )
        try:
            run = Run.load(self.workspace.run_folder(beam_run))
        except FileNotFoundError:
            raise FileNotFoundError(
                f"no finished run {beam_run!r} in {self.workspace.runs_folder}"
            ) from None
        if run.run_type is not RunType.RUN_TO_SAMPLE:
            raise ValueError(f"{beam_run!r} is a {run.run_type} run, not a run_to_sample")
        if folder_name == beam_run:
            raise ValueError(f"folder_name {folder_name!r} is the run to the sample it starts from")
        return run

    def get_scan_points(
        self, variable: str, values: ArrayLike, instr_params: Mapping[str, float]
    ) -> list[dict[str, float]]:
        """Return the points of a scan.

        Parameters
        ----------
        variable:
            Instrument parameter to scan.
        values:
            Values of ``variable``.
        instr_params:
            The other instrument parameters, the same for every point.

        Returns
        -------
        :
            ``{**instr_params, variable: value}`` for each value.

        Raises
        ------
        KeyError
            If a name is not an instrument parameter, before any run.
        """
        names = {p.name for p in self.instrument.instr_params}
        if unknown := sorted({variable, *instr_params} - names):
            raise KeyError(f"no instrument parameter {unknown}, only {sorted(names)}")
        return [
            {**instr_params, variable: value} for value in np.asarray(values, dtype=float).tolist()
        ]

    def run(
        self,
        folder_name: str | None = None,
        ncount: float = 1e7,
        seed: int | None = None,
        mpi: int | str = "auto",
        custom_flags: str | None = None,
        **instr_params: float,
    ) -> Run:
        """Simulate source -> detector once.

        Parameters
        ----------
        folder_name:
            Run folder name. If ``None``, ``<time>_run``.
        ncount:
            Number of neutrons.
        seed:
            Random seed. If ``None``, McStas picks one.
        mpi:
            Number of MPI processes, or ``"auto"``.
        custom_flags:
            Extra flags for mcrun.
        **instr_params:
            Instrument parameters by McStas name.

        Returns
        -------
        :
            The run, with its file ``data_file``, ``<runs folder>/<folder>/mccode.h5``.
        """
        run, inst = self.prepare(
            folder_name,
            RunType.RUN,
            instr_params,
            ncount=ncount,
            seed=seed,
            mpi=mpi,
            custom_flags=custom_flags,
        )
        return self._run_once(run, inst)

    def scan(
        self,
        variable: str,
        values: ArrayLike,
        folder_name: str | None = None,
        ncount: float = 1e7,
        n_jobs: int | None = None,
        seed: int | None = None,
        custom_flags: str | None = None,
        **instr_params: float,
    ) -> Run:
        """Simulate source -> detector for every value of an instrument parameter.

        The values run in parallel as single-process mcrun jobs.

        Parameters
        ----------
        variable:
            Instrument parameter to scan, e.g. ``"b_rot"``.
        values:
            Values of ``variable``, e.g. ``np.linspace(-0.5, -0.3, 11)``.
        folder_name:
            Run folder name. If ``None``, ``<time>_scan_<variable>``.
        ncount:
            Number of neutrons per value.
        n_jobs:
            Number of jobs at a time. If ``None``, the number of cores.
        seed:
            Random seed; value ``n`` gets ``seed + n + 1``, so the jobs do not repeat each
            other's random numbers. If ``None``, a random seed.
        custom_flags:
            Extra flags for mcrun.
        **instr_params:
            The other instrument parameters by McStas name.

        Returns
        -------
        :
            The run, with its file ``data_file``, ``<runs folder>/<folder>/mccode.h5``, which
            has one entry per value.

        Raises
        ------
        KeyError
            If a name is not an instrument parameter, before anything runs.
        """
        points = self.get_scan_points(variable, values, instr_params)  # checked before anything
        run, inst = self.prepare(
            folder_name,
            RunType.SCAN,
            instr_params,
            ncount=ncount,
            seed=random_seed() if seed is None else seed,  # picked here, to record it
            variable=variable,
            values=[point[variable] for point in points],
            mpi=None,  # no MPI: the points run in parallel instead
            custom_flags=custom_flags,
        )
        return self._run_points(run, inst, points, n_jobs)

    def run_to_sample(
        self,
        folder_name: str | None = None,
        ncount: float = 1e8,
        seed: int | None = None,
        mpi: int | str = "auto",
        **instr_params: float,
    ) -> Run:
        """Simulate source -> ``sample_pos`` and dump the beam to an MCPL file.

        This is part 1 of a split run: run it once, with a large ``ncount``. The run folder also
        holds a ``mccode.h5`` with the incident monitor.

        The beam file only fits instruments with the same components up to ``sample_pos``, and the
        same instrument parameters they use. Make a new one after changing any of them.

        Parameters
        ----------
        folder_name:
            Run folder name. If ``None``, ``<time>_run_to_sample``.
        ncount:
            Number of neutrons.
        seed:
            Random seed. If ``None``, McStas picks one.
        mpi:
            Number of MPI processes, or ``"auto"``.
        **instr_params:
            Instrument parameters by McStas name.

        Returns
        -------
        :
            The run, with its beam ``beam_file``, ``<runs folder>/<folder>/beam.mcpl.gz``. It
            also becomes :attr:`beam_run`.

        Raises
        ------
        ValueError
            If the instrument has no ``sample_pos``, before anything runs.
        """
        self._check_sample_pos()  # before anything
        run, inst = self.prepare(
            folder_name, RunType.RUN_TO_SAMPLE, instr_params, ncount=ncount, seed=seed, mpi=mpi
        )
        inst.run_to(
            SAMPLE_POS,
            run_name="beam",
            comment=f"{ncount:g} neutrons",
            filename='"beam.mcpl"',
            # the ray count is stored in the file and the weights scaled by it; McStas 3.9 would
            # default to the raw weights of mode 2, which are deprecated
            weight_mode=1,
        )
        run = self._run_once(run, inst)
        self.beam_run = run.name
        return run

    def run_from_sample(
        self,
        folder_name: str | None = None,
        beam_run: str | None = None,
        seed: int | None = None,
        mpi: int | str = "auto",
        **instr_params: float,
    ) -> Run:
        """Simulate ``sample_pos`` -> sample -> detector once, from a beam file.

        This is part 2 of a split run, and may run any time later, even in another session. The
        number of neutrons is the number in the beam file.

        Parameters
        ----------
        folder_name:
            Run folder name. If ``None``, ``<time>_run_from_sample``.
        beam_run:
            Name of the run to the sample whose beam to start from, in the runs folder of the
            workspace. If ``None``, :attr:`beam_run`, the last one of this simulation.
        seed:
            Random seed. If ``None``, McStas picks one.
        mpi:
            Number of MPI processes, or ``"auto"``. Each reads its own part of the beam file.
        **instr_params:
            Instrument parameters by McStas name. Those the beam depends on must be the ones of the
            beam file.

        Returns
        -------
        :
            The run, with its file ``data_file``, ``<runs folder>/<folder>/mccode.h5``.

        Raises
        ------
        ValueError
            If the instrument has no ``sample_pos``, before anything runs.
        FileNotFoundError
            If there is no run to the sample to start from, before anything runs.
        """
        run, inst = self._prepare_from_sample(
            folder_name, RunType.RUN_FROM_SAMPLE, instr_params, beam_run, seed=seed, mpi=mpi
        )
        return self._run_once(run, inst)

    def scan_from_sample(
        self,
        variable: str,
        values: ArrayLike,
        folder_name: str | None = None,
        beam_run: str | None = None,
        n_jobs: int | None = None,
        seed: int | None = None,
        **instr_params: float,
    ) -> Run:
        """Simulate ``sample_pos`` -> detector from a beam file, for every value of a parameter.

        This is part 2 of a split run, for a scan. The values run in parallel as single-process
        mcrun jobs, all starting from the same beam file.

        Parameters
        ----------
        variable:
            Instrument parameter to scan, e.g. ``"b_rot"``.
        values:
            Values of ``variable``, e.g. ``np.linspace(-0.5, -0.3, 11)``.
        folder_name:
            Run folder name. If ``None``, ``<time>_scan_<variable>_from_sample``.
        beam_run:
            Name of the run to the sample whose beam to start from, in the runs folder of the
            workspace. If ``None``, :attr:`beam_run`, the last one of this simulation.
        n_jobs:
            Number of jobs at a time. If ``None``, the number of cores.
        seed:
            Random seed; value ``n`` gets ``seed + n + 1``, so the jobs do not repeat each
            other's random numbers. If ``None``, a random seed.
        **instr_params:
            The other instrument parameters by McStas name.

        Returns
        -------
        :
            The run, with its file ``data_file``, ``<runs folder>/<folder>/mccode.h5``, which
            has one entry per value.

        Raises
        ------
        KeyError
            If a name is not an instrument parameter, before anything runs.
        ValueError
            If the instrument has no ``sample_pos``, before anything runs.
        FileNotFoundError
            If there is no run to the sample to start from, before anything runs.
        """
        points = self.get_scan_points(variable, values, instr_params)  # checked before anything
        run, inst = self._prepare_from_sample(
            folder_name,
            RunType.SCAN_FROM_SAMPLE,
            instr_params,
            beam_run,
            seed=random_seed() if seed is None else seed,  # picked here, to record it
            variable=variable,
            values=[point[variable] for point in points],
            mpi=None,  # no MPI: the points run in parallel instead
        )
        return self._run_points(run, inst, points, n_jobs)
