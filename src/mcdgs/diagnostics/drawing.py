"""Describe and draw any instrument: its printout, its McStas instrument file and its views.

- :func:`write_printout` writes ``print(instrument)``: its parameters and the position of every
  component.
- :func:`write_instr_file` writes the McStas instrument file into its build folder, and copies
  it to where the other files go.
- :func:`draw_instrument` draws the instrument with McStasScript's ``show_instrument``, which runs
  mcdisplay, in any projection, and :func:`component_range` gives the components a view draws.
- :func:`draw_diagram` draws McStasScript's diagram of how the components are placed, and
  :func:`draw_intensity_diagram` the same with the intensity before each component, from a run.
- :func:`remove_mcdisplay_folders` clears what mcdisplay leaves in the build folder.

McStasScript's drawings end in ``plt.show()`` and return None, so they are drawn with
``plt.show`` held back, and the figures they would show are saved instead.
"""

import shutil
from pathlib import Path
from typing import Any
from unittest.mock import patch

import matplotlib.pyplot as plt
from mcstasscript.interface.instr import McStas_instr

from mcdgs.core import Instrument
from mcdgs.diagnostics import cli


def write_printout(instrument: Instrument, path: Path) -> Path:
    """Write ``print(instrument)``: its parameters and the position of every component.

    Parameters
    ----------
    instrument:
        The instrument.
    path:
        The file, such as ``trex_print.txt``.

    Returns
    -------
    :
        ``path``.
    """
    return cli.save_text(f"{instrument}\n", path)


def write_instr_file(mcstas_instr: McStas_instr, folder: Path) -> Path:
    """Write the McStas instrument file into the build folder, and copy it to ``folder``.

    Parameters
    ----------
    mcstas_instr:
        The McStasScript instrument, as :func:`.build_mcstas` builds it.
    folder:
        The folder of the copy.

    Returns
    -------
    :
        The copy, ``folder/<instrument name>.instr``.
    """
    mcstas_instr.write_full_instrument()  # into its input_path, the build folder
    written = Path(mcstas_instr.input_path) / f"{mcstas_instr.name}.instr"
    copy = Path(shutil.copyfile(written, folder / written.name))
    print(f"Saved {copy}")
    return copy


def component_range(instrument: Instrument, after: str, before: str) -> tuple[int, int]:
    """Return the indices of the components between two, as ``show_instrument`` takes them.

    A view along the beam leaves out such components as a wide source and a detector around
    the sample, which would cover the rest.

    Parameters
    ----------
    instrument:
        The instrument.
    after:
        The component just before the first one drawn.
    before:
        The component just after the last one drawn.

    Returns
    -------
    :
        The index of the first component after ``after``, and that of ``before``, the first
        one not drawn: ``index_min`` and ``index_max`` of ``show_instrument``.

    Raises
    ------
    ValueError
        If the instrument has no component of either name.
    """
    names = list(instrument.components)
    if missing := [name for name in (after, before) if name not in names]:
        raise ValueError(f"{instrument.name} has no component {missing}")
    return names.index(after) + 1, names.index(before)


def draw_instrument(
    mcstas_instr: McStas_instr,
    path: Path,
    *,
    dpi: int = 150,
    equal_aspect: bool = True,
    **view: Any,
) -> Path:
    """Draw the instrument with McStasScript's ``show_instrument``, and save the figure.

    mcdisplay draws every component, whatever its ``WHEN``. The folder it leaves in the build
    folder is removed once the view is drawn.

    Parameters
    ----------
    mcstas_instr:
        The McStasScript instrument.
    path:
        The file of the figure.
    dpi:
        Its resolution, in dots per inch.
    equal_aspect:
        True gives both axes the same scale, which ``show_instrument`` has no option for.
    **view:
        The view, as ``show_instrument`` takes it, in its ``matplotlib_2d`` backend: such as
        ``projection="zx"``, from above with z along the beam, or ``"xy"``, along the beam;
        ``limits`` in m, ``[min, max]`` across then up; ``index_min`` and ``index_max``
        (:func:`component_range`); and ``width`` and ``height``, in pixels at 100 dpi.

    Returns
    -------
    :
        ``path``.
    """
    with patch.object(plt, "show"):
        mcstas_instr.show_instrument(backend="matplotlib_2d", **view)
    remove_mcdisplay_folders(mcstas_instr.input_path)  # not needed once drawn
    if equal_aspect:
        plt.gca().set_aspect("equal")
    return cli.save_figure(plt.gcf(), path, dpi=dpi)


def draw_diagram(mcstas_instr: McStas_instr, legend_path: Path, diagram_path: Path) -> list[Path]:
    """Draw McStasScript's diagram of how the components are placed and its legend, and save them.

    ``show_diagram()`` lists the components in beam order, each with an arrow to the component
    it is placed ``AT ... RELATIVE`` to, in blue, and ``ROTATED ... RELATIVE`` to, in red. It
    draws its legend first, in a figure of its own: the color of each kind of McStas component,
    and the dashed box of a component with a ``WHEN``. Both are saved without their empty
    margins.

    Parameters
    ----------
    mcstas_instr:
        The McStasScript instrument.
    legend_path:
        The file of the legend.
    diagram_path:
        The file of the diagram.

    Returns
    -------
    :
        ``legend_path`` and ``diagram_path``.
    """
    plt.close("all")  # the two figures it draws are then the only ones
    with patch.object(plt, "show"):
        mcstas_instr.show_diagram()
    legend, diagram = (plt.figure(number) for number in plt.get_fignums())  # in drawing order
    paths = [legend_path, diagram_path]
    for figure, path in zip((legend, diagram), paths, strict=True):
        figure.savefig(path, dpi=100, bbox_inches="tight")
        print(f"Saved {path}")
    plt.close("all")
    return paths


def draw_intensity_diagram(
    mcstas_instr: McStas_instr,
    path: Path,
    intensity_limits: tuple[float, float] | None = None,
) -> Path:
    """Run the instrument and draw McStasScript's diagram with the intensity before each component.

    ``show_diagram(analysis=True)`` runs McStas with the settings of ``mcstas_instr``, into its
    ``output_path``, with a monitor before every component. It then draws the diagram of
    :func:`draw_diagram`, and next to it the intensity just before each component, in n/s, and
    the number of rays, on a log scale.

    Parameters
    ----------
    mcstas_instr:
        The McStasScript instrument, with the settings of the run.
    path:
        The file of the figure.
    intensity_limits:
        The range of the intensity axis, in n/s, which ``show_diagram`` has no option for: its
        ``limits`` are those of a ``variable``. The ray count keeps its ratio to the
        intensity, so both share their ticks. If ``None``, as McStasScript draws it.

    Returns
    -------
    :
        ``path``.
    """
    with patch.object(plt, "show"):
        mcstas_instr.show_diagram(analysis=True)
    diagram = plt.gcf()
    if intensity_limits is not None:  # the ray count zoomed with it, keeping its ratio
        intensity_ax = next(ax for ax in diagram.axes if ax.get_xlabel() == "Intensity [n/s]")
        rays_ax = next(ax for ax in diagram.axes if ax.get_xlabel() == "Ray count")
        rays_per_intensity = rays_ax.get_xlim()[1] / intensity_ax.get_xlim()[1]
        low, high = intensity_limits
        intensity_ax.set_xlim(low, high)
        rays_ax.set_xlim(low * rays_per_intensity, high * rays_per_intensity)
    return cli.save_figure(diagram, path, dpi=100)


def remove_mcdisplay_folders(build_folder: str | Path) -> list[Path]:
    """Remove the folders mcdisplay left in the build folder of an instrument.

    McStasScript's ``show_instrument`` runs mcdisplay, which writes a new folder
    ``<instrument>_mcdisplay``, ``<instrument>_mcdisplay_0``, ... beside the instrument file at
    every call. They are not needed once the instrument is drawn.

    Parameters
    ----------
    build_folder:
        The build folder of the instrument, McStasScript's ``input_path``.

    Returns
    -------
    :
        The removed folders.
    """
    folders = sorted(p for p in Path(build_folder).glob("*_mcdisplay*") if p.is_dir())
    for folder in folders:
        shutil.rmtree(folder)
    return folders
