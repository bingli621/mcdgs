"""Run McStasScript's diagnostics of the beam, with the folders and records of a :class:`.Simulation`.

McStasScript's ``Diagnostics`` records the rays at points along the instrument and draws them in
views: histograms of one ray variable, or of one against another. It runs McStas itself, so a
:class:`McStasScriptDiagnostics` keeps its points and views, and on every ``run`` gives them to
a new ``Diagnostics`` of the run that :meth:`.Simulation.prepare` makes.
"""

from typing import Any, Literal

import mcstasscript as ms

from mcdgs.core import Instrument
from mcdgs.mcstas.runs import Run, RunType
from mcdgs.mcstas.simulation import Simulation
from mcdgs.mcstas.workspace import Workspace


class McStasScriptDiagnostics:
    """McStasScript's diagnostics of the beam of an instrument.

    Add the points and the views, then ``run``, as often as needed, and ``plot`` the last run.
    Every run builds the instrument anew, as a run of :class:`.Simulation` does, so a change to
    the instrument applies to the next run.

    Parameters
    ----------
    instrument:
        The instrument to diagnose.
    workspace:
        The folders to write into, as for :class:`.Simulation`.

    Attributes
    ----------
    diagnostics:
        McStasScript's ``Diagnostics`` of the last run, with its rays. ``None`` before the first
        run.
    """

    def __init__(self, instrument: Instrument, workspace: Workspace | None = None) -> None:
        self.simulation = Simulation(instrument, workspace)
        self.diagnostics: ms.Diagnostics | None = None
        self._points: list[dict[str, Any]] = []
        self._views: list[tuple[str, str | None, dict[str, Any]]] = []

    def add_point(
        self,
        before: str | None = None,
        after: str | None = None,
        rays: int | Literal["all"] = "all",
    ) -> None:
        """Record the rays just before or just after a component.

        Parameters
        ----------
        before, after:
            The name of the component, one of the two or both.
        rays:
            Number of rays to record. ``"all"``, as McStasScript's limited list stops the run
            early, which leaves the intensities too low.

        Raises
        ------
        ValueError
            If neither ``before`` nor ``after`` is given.
        KeyError
            If a component is not one of the instrument.
        """
        if before is None and after is None:
            raise ValueError("give the component the point is before or after")
        names = self.simulation.instrument.components
        given = {name for name in (before, after) if name is not None}
        if unknown := sorted(given - set(names)):
            raise KeyError(f"no component {unknown} in {self.simulation.instrument.name}")
        self._points.append({"before": before, "after": after, "rays": rays})

    def add_view(self, axis1: str, axis2: str | None = None, **kwargs: Any) -> None:
        """Draw every point as a histogram of one ray variable, or of one against another.

        Parameters
        ----------
        axis1, axis2:
            The ray variables, e.g. ``"x"``, ``"y"`` or ``"l"``, the wavelength. A histogram of
            ``axis1`` alone if ``axis2`` is ``None``.
        **kwargs:
            The other options of McStasScript's ``add_view``, e.g. ``bins``. A view of two
            variables is drawn with ``colormap="turbo"`` unless another is given.
        """
        defaults = {} if axis2 is None else {"colormap": "turbo"}
        self._views.append((axis1, axis2, defaults | kwargs))

    def run(
        self,
        folder_name: str | None = None,
        ncount: float = 1e6,
        seed: int | None = None,
        mpi: int | str = "auto",
        **instr_params: float,
    ) -> Run:
        """Simulate source -> detector once, recording the rays at the points.

        The run writes McStas's own format, which McStasScript reads: the record names a
        ``data_file`` that does not exist.

        Parameters
        ----------
        folder_name:
            Run folder name. If ``None``, ``<time>_diagnose``.
        ncount:
            Number of neutrons. The points keep their rays, so fewer than for
            :meth:`.Simulation.run`.
        seed:
            Random seed. If ``None``, McStas picks one.
        mpi:
            Number of MPI processes, or ``"auto"``.
        **instr_params:
            Instrument parameters by McStas name.

        Returns
        -------
        :
            The run. Its rays are in :attr:`diagnostics`.

        Raises
        ------
        ValueError
            If there is no point or no view, before anything runs.
        """
        if not self._points or not self._views:
            raise ValueError("add a point and a view before the run")
        run, inst = self.simulation.prepare(
            folder_name,
            RunType.DIAGNOSE,
            instr_params,
            ncount=ncount,
            seed=seed,
            nexus=False,
            mpi=mpi,
        )
        diagnostics = ms.Diagnostics(inst)  # runs a copy of inst, with its settings
        for point in self._points:
            diagnostics.add_point(**point)
        for axis1, axis2, kwargs in self._views:
            diagnostics.add_view(axis1, axis2, **kwargs)
        diagnostics.run()  # compiles in its build folder
        self.diagnostics = diagnostics
        return run.finish()

    def plot(self) -> None:
        """Draw the last run, a row per point and a column per view, on a new figure.

        McStasScript does not return the figure: ``plt.gcf()`` gives it.

        Raises
        ------
        RuntimeError
            If there has been no run.
        """
        if self.diagnostics is None:
            raise RuntimeError("run the diagnostics before plotting them")
        self.diagnostics.plot()
