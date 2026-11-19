"""Run McStasScript's diagnostics of the beam, with the folders and records of a :class:`.Simulation`.

McStasScript's ``Diagnostics`` records the rays at points along the instrument and draws them in
views: histograms of one ray variable, or of one against another. Its ``show_diagram`` with
``analysis=True`` draws the components with the intensity before each. Both run McStas
themselves, so a :class:`McStasScriptDiagnostics` hands them the run that
the :class:`.Simulation` prepares, and finishes its record once they have run. It keeps the
points and views, and on every ``run`` gives them to a new ``Diagnostics``.
"""

from typing import Any, Literal
from unittest.mock import patch

import matplotlib.pyplot as plt
import mcstasscript as ms
from matplotlib.figure import Figure

from mcdgs.core import Instrument
from mcdgs.mcstas.runs import Run, RunType
from mcdgs.mcstas.simulation import Simulation
from mcdgs.mcstas.workspace import Workspace


class McStasScriptDiagnostics:
    """McStasScript's diagnostics of the beam of an instrument.

    Add the points and the views, then ``run``, as often as needed, and ``plot`` the last run.
    :meth:`run_analysis` runs and draws the diagram of the intensity, without them. Every
    run builds the instrument anew, as a run of :class:`.Simulation` does, so a change to
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

    @property
    def workspace(self) -> Workspace:
        """The folders the runs write into: those of :attr:`simulation`."""
        return self.simulation.workspace

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
        run, inst = self.simulation._prepare(
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

    def plot(self) -> Figure:
        """Draw the last run, a row per point and a column per view, on a new figure.

        Returns
        -------
        :
            The figure.

        Raises
        ------
        RuntimeError
            If there has been no run.
        """
        if self.diagnostics is None:
            raise RuntimeError("run the diagnostics before plotting them")
        self.diagnostics.plot()  # draws on a new figure, which it does not return
        return plt.gcf()

    def run_analysis(
        self,
        folder_name: str | None = None,
        ncount: float = 1e7,
        seed: int | None = None,
        mpi: int | str = "auto",
        intensity_limits: tuple[float, float] | None = None,
        **instr_params: float,
    ) -> tuple[Run, Figure]:
        """Simulate once, and draw McStasScript's diagram with the intensity before each component.

        ``show_diagram(analysis=True)`` runs McStas with a monitor before every component, then
        draws the diagram of the components, with arrows from each to the one it is placed
        relative to, and next to it the intensity just before each component, in n/s, and the
        number of rays, on a log scale. The points and views play no part. The run writes
        McStas's own format, as :meth:`run` does, and its folder holds the monitors of the
        instrument too.

        Parameters
        ----------
        folder_name:
            Run folder name. If ``None``, ``<time>_analysis``.
        ncount:
            Number of neutrons.
        seed:
            Random seed. If ``None``, McStas picks one.
        mpi:
            Number of MPI processes, or ``"auto"``.
        intensity_limits:
            The range of the intensity axis, in n/s, which ``show_diagram`` has no option for:
            its ``limits`` are those of a ``variable``. The ray count keeps its ratio to the
            intensity, so both share their ticks. If ``None``, as McStasScript draws it.
        **instr_params:
            Instrument parameters by McStas name.

        Returns
        -------
        :
            The run, and the figure of the diagram.
        """
        run, inst = self.simulation._prepare(
            folder_name,
            RunType.ANALYSIS,
            instr_params,
            ncount=ncount,
            seed=seed,
            nexus=False,
            mpi=mpi,
        )
        with patch.object(plt, "show"):  # McStasScript shows the figure and returns None
            inst.show_diagram(analysis=True)
        diagram = plt.gcf()
        if intensity_limits is not None:  # the ray count zoomed with it, keeping its ratio
            intensity_ax = next(ax for ax in diagram.axes if ax.get_xlabel() == "Intensity [n/s]")
            rays_ax = next(ax for ax in diagram.axes if ax.get_xlabel() == "Ray count")
            rays_per_intensity = rays_ax.get_xlim()[1] / intensity_ax.get_xlim()[1]
            low, high = intensity_limits
            intensity_ax.set_xlim(low, high)
            rays_ax.set_xlim(low * rays_per_intensity, high * rays_per_intensity)
        return run.finish(), diagram
