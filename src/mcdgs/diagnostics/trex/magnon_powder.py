"""The S(Q, E) of a powder magnon in every sub-pulse of T-Rex: the resolution across the RRM.

:class:`TRexMagnonPowder` is :class:`.TRex` with the powder of :func:`.magnon_powder` on
``sample_pos``, with the beamstop and the banana detector, as
:mod:`.diagnostics.trex.elastic_resolution_incoherent` has its vanadium: the same cylinder, the
same choppers, so the same 15 sub-pulses, each with its own incident energy. Where vanadium
gives one line at E = 0, the powder gives a dispersive magnon, from its gap of
:data:`.MAGNON_GAP` at small |Q| to the top of its band, :data:`.MAGNON_TOP`. The S(|Q|, E) of
the powder, from :func:`.write_magnon_sqw`, is a band at each |Q| between sharp edges, with a line
far narrower than the resolution: how much each sub-pulse blurs those edges is its resolution,
across E and |Q|. The powder scatters once, within ``d_phi`` of the horizontal plane: the height
of the window of :func:`.solid_angles.detector_window`, which covers the banana detector.

It runs in the mode of ``MODES``, ``"FF"``, the wide slits of the choppers, for high flux; a
caller may pass others. :func:`simulate` writes the S(|Q|, E) file of the powder into the builds
folder of the workspace, :func:`sqw_path`, then runs each mode in two parts, with
:class:`.Simulation`, as the vanadium does: from the source to ``sample_pos``, then from the beam
there through the powder to the detector, into the folders of :func:`run_names`.

The reduction is that of the vanadium (:func:`.elastic_resolution_incoherent.reduce_runs`): the
sub-pulses, the pixels past ``MIN_TWO_THETA`` from the beam, the neutrons of each sub-pulse
normalized by the solid angle of their pixel and by its flux, and :func:`.q_e_maps`, S(Q, E)
normalized by the coverage of each bin. :func:`plot_magnon_maps` draws the 15 maps in one figure,
with the edges of the band of the model, :func:`.magnon_band`, over them.

:func:`main` is :func:`simulate`, which runs McStas, then :func:`analyze`, which reduces and
draws. Run it to do both::

    pixi run python -c "from mcdgs.diagnostics.trex.magnon_powder import main; main()"

Running the module only reduces the runs again, without McStas, with :func:`analyze`::

    pixi run python -m mcdgs.diagnostics.trex.magnon_powder
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipp as sc
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from mcdgs.core.components import Isotropic_Sqw
from mcdgs.diagnostics import cli
from mcdgs.diagnostics import elastic_resolution as resolution
from mcdgs.diagnostics.trex import elastic_resolution_incoherent as vanadium
from mcdgs.diagnostics.trex import sample_runs as runs
from mcdgs.diagnostics.trex.solid_angles import detector_window
from mcdgs.instruments.trex.choppers import Mode
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.mcstas import Simulation, Workspace
from mcdgs.reduction.inelastic import normalized_events
from mcdgs.samples import magnon_powder, write_magnon_sqw
from mcdgs.samples.inelastic.powder import MAGNON_TOP, magnon_band

NAME = "trex_magnon_powder"
SQW_NAME = "magnon_powder.sqw"  # in the builds folder of the workspace
MODES: tuple[Mode, ...] = ("FF",)  # high flux
P_INTERACT = 0.9  # the fraction of the neutrons made to scatter in the powder
SAMPLE_SPLIT = 4  # each neutron of the beam goes through the powder this many times
NCOUNT = 5e9  # McStas rays of the first part, from the masked source
ENERGY_STEP = 0.25  # meV, the bin width in E of the maps
Q_STEP = 0.1  # 1/Å, the bin width in |Q| of the maps
# the E of the maps: the magnon is created only, so from a little below 0 to past the band
ENERGY_RANGE = (-2.0, MAGNON_TOP + 3.0)  # meV
EDGE_COLOR = "white"  # of the edges of the band of the model, over the maps


def sqw_path(workspace: Workspace | None = None) -> Path:
    """Return the S(|Q|, E) file of the powder: ``SQW_NAME`` in the builds folder.

    McStas reads it from the folder it runs in, so the path is absolute.

    Parameters
    ----------
    workspace:
        The workspace. If ``None``, that of :meth:`.Workspace.find`.
    """
    workspace = Workspace.find() if workspace is None else workspace
    return (workspace.builds_folder / SQW_NAME).resolve()


# --8<-- [start:instrument]
def sample(sqw_file: str | Path | None = None) -> Isotropic_Sqw:
    """Return the powder of :func:`.magnon_powder`, scattering once towards the detector.

    It scatters within ``d_phi`` of the horizontal plane, the height of the window of
    :func:`.solid_angles.detector_window`, and makes ``P_INTERACT`` of the neutrons scatter. Each
    neutron goes through it ``SAMPLE_SPLIT`` times (McStas SPLIT), to fill S(Q, E) more smoothly.

    Parameters
    ----------
    sqw_file:
        Its S(|Q|, E) file. If ``None``, that of :func:`sqw_path`.
    """
    return magnon_powder(
        sqw_path() if sqw_file is None else sqw_file,
        name="sample",
        order=1,
        d_phi=detector_window()["focus_ah"],
        p_interact=P_INTERACT,
        split=SAMPLE_SPLIT,
    )


class TRexMagnonPowder(TRex):
    """T-Rex with a powder of a magnon, from the source to the banana detector.

    The powder is that of :func:`sample`.

    Parameters
    ----------
    params:
        Control parameters.
    source_masked:
        As in :class:`.TRex`. True by default: the masked source of
        :func:`.sources.masked_source` spends the rays where the choppers let neutrons through.
    sqw_file:
        As :func:`sample`.
    """

    name = "trex_magnon_powder"  # the McStas instrument name

    def __init__(
        self,
        params: Params = DEFAULT_PARAMS,
        source_masked: bool = True,
        sqw_file: str | Path | None = None,
    ) -> None:
        super().__init__(sample=sample(sqw_file), params=params, source_masked=source_masked)


def define_instrument(mode: Mode = "FF") -> TRexMagnonPowder:
    """Return the instrument of this diagnostic: :class:`TRexMagnonPowder` at its defaults.

    Parameters
    ----------
    mode:
        The slits of the pulse-shaping and monochromatic choppers: ``"FF"``, the default, for
        high flux, ``"RR"`` for high resolution, or ``"RF"`` or ``"FR"`` between.
    """
    params = Params(mode=mode)  # the default control parameters, in the mode
    return TRexMagnonPowder(params=params)


# --8<-- [end:instrument]


def run_names(mode: Mode) -> tuple[str, str]:
    """Return the run folders of a mode: the part to ``sample_pos``, and the part from it.

    Examples
    --------
    >>> run_names("RR")[1]
    'trex_magnon_powder_rr_run_from_sample'
    """
    return f"{NAME}_{mode.lower()}_run_to_sample", f"{NAME}_{mode.lower()}_run_from_sample"


def simulate(
    ncount: float = NCOUNT, seed: int = 1, modes: tuple[Mode, ...] = MODES
) -> dict[Mode, tuple[Path, Path]]:
    """Write the S(|Q|, E) of the powder, then run :class:`TRexMagnonPowder` in each mode.

    The file is that of :func:`.write_magnon_sqw` at its defaults, at :func:`sqw_path`: McStas
    reads it when the instrument starts, in both parts. In each mode, the first part runs from
    the source to ``sample_pos``, the second from that beam through the powder to the detector,
    both into the folders of :func:`run_names`. Each replaces the folder of an earlier run.

    Parameters
    ----------
    ncount:
        Number of McStas rays of the first part, from the masked source. The second part runs
        every neutron of the beam file.
    seed:
        Seed of the random numbers of the first part, in every mode; the second part takes
        ``seed + 1``.
    modes:
        The modes to run.

    Returns
    -------
    :
        The folders of the first and the second part, by mode.
    """
    write_magnon_sqw(sqw_path())
    folders = {}
    for mode in modes:
        to_name, from_name = run_names(mode)
        sim = Simulation(define_instrument(mode))
        to_sample = sim.run_to_sample(folder_name=to_name, ncount=ncount, seed=seed)
        from_sample = sim.run_from_sample(folder_name=from_name, seed=seed + 1)
        folders[mode] = (to_sample.folder, from_sample.folder)
    return folders


def reduce_runs(
    mode: Mode = "FF", to_sample: str | Path | None = None, from_sample: str | Path | None = None
) -> tuple[sc.DataArray, sc.DataArray]:
    """Reduce the two parts of a McStas run of :func:`simulate` to S(Q, E), without McStas.

    As :func:`.elastic_resolution_incoherent.reduce_runs`, with the bins ``ENERGY_STEP`` and
    ``Q_STEP``.

    Parameters
    ----------
    mode:
        The mode of the run.
    to_sample, from_sample:
        The folders of the first and the second part. If ``None``, those of :func:`run_names`
        in the runs folder of :meth:`.Workspace.find`, where :func:`simulate` writes.

    Returns
    -------
    :
        The sub-pulses of :func:`.sample_runs.subpulse_timing`, and the maps of
        :func:`.q_e_maps`, from the pixels of :func:`.scattered_pixels`.
    """
    workspace = Workspace.find()
    to_name, from_name = run_names(mode)
    to_sample = Path(to_sample or workspace.run_folder(to_name))
    from_sample = Path(from_sample or workspace.run_folder(from_name))
    pulses = runs.subpulse_timing_of_run(to_sample, define_instrument(mode))
    binned = resolution.scattered_pixels(runs.read_detector(from_sample), vanadium.MIN_TWO_THETA)
    pixels = runs.pixel_solid_angle(binned)
    events = normalized_events(binned, pulses, pixels)
    maps = resolution.q_e_maps(
        events,
        pulses,
        resolution.two_theta(binned),
        pixels,
        energy_step=ENERGY_STEP,
        q_step=Q_STEP,
    )
    return pulses, maps


def crop(maps: sc.DataArray, energy_range: tuple[float, float] = ENERGY_RANGE) -> sc.DataArray:
    """Keep the bins of E within ``energy_range``, in meV, of the maps of :func:`.q_e_maps`.

    Examples
    --------
    >>> edges = sc.linspace("energy_transfer", -10.0, 10.0, 21, unit="meV")
    >>> maps = sc.DataArray(
    ...     sc.zeros(dims=["energy_transfer"], shape=[20]), coords={"energy_transfer": edges}
    ... )
    >>> cropped = crop(maps, (-2.0, 3.0))
    >>> cropped.coords["energy_transfer"].values[[0, -1]].tolist()
    [-2.0, 3.0]
    """
    low, high = (sc.scalar(value, unit="meV") for value in energy_range)
    unit = maps.coords["energy_transfer"].unit
    return maps["energy_transfer", low.to(unit=unit) : high.to(unit=unit)]


def plot_magnon_maps(maps: sc.DataArray, title: str) -> Figure:
    """Draw the S(Q, E) of every sub-pulse, with the edges of the band of the model over it.

    The maps are those of :func:`.plot_q_e_maps`, within ``ENERGY_RANGE``: one panel per
    sub-pulse, from the highest Ei, on one color scale, white where the sub-pulse cannot reach.
    Over each, the lowest and the highest energy of the magnon at each |Q|, :func:`.magnon_band`,
    are dashed lines: the band of the powder has sharp edges, and how far the intensity spreads
    past them is the resolution of the sub-pulse.

    Parameters
    ----------
    maps:
        The maps of :func:`.q_e_maps`.
    title:
        The title of the figure.

    Returns
    -------
    :
        The figure.
    """
    maps = crop(maps)
    fig = resolution.plot_q_e_maps(maps, title)
    q_edges = maps.coords["Q"].to(unit="1/Å").values
    q = np.linspace(q_edges[0], q_edges[-1], 400)
    low, high = magnon_band(q)
    panels = [ax for ax in fig.axes[: maps.sizes["subpulse"]] if ax.get_visible()]
    for ax in panels:
        for edge in (low, high):
            ax.plot(q, edge, color=EDGE_COLOR, lw=0.8, ls="--")
    if panels:
        panels[0].legend(
            [Line2D([], [], color=EDGE_COLOR, lw=0.8, ls="--")],
            ["Edges of the band, without resolution"],
            loc="upper right",
            fontsize="small",
            facecolor="0.3",
            labelcolor=EDGE_COLOR,
        )
    return fig


def analyze(out: Path | None = None, modes: tuple[Mode, ...] = MODES) -> list[Path]:
    """Reduce the runs of each mode and draw the S(Q, E) of all its sub-pulses.

    It runs no McStas, so it can run again on the same runs; :func:`main` runs them first. For
    each mode, ``m`` its name in lower case, the file is ``trex_magnon_powder_m_sqe.png``, from
    :func:`plot_magnon_maps`.

    Parameters
    ----------
    out:
        Folder of the files. If ``None``, ``trex_magnon_powder`` in the results folder of
        :meth:`.Workspace.find`.
    modes:
        The modes to reduce, whose runs :func:`simulate` made.

    Returns
    -------
    :
        The paths of the files.
    """
    folder = cli.output_folder(out, NAME)
    files = []
    for mode in modes:
        pulses, maps = reduce_runs(mode)
        conditions = define_instrument(mode).running_conditions()
        fig = plot_magnon_maps(
            maps,
            f"T-Rex, a powder magnon: S(Q, E) of each of the {pulses.sizes['subpulse']} "
            f"sub-pulses\n{conditions}",
        )
        files.append(cli.save_figure(fig, folder / f"{NAME}_{mode.lower()}_sqe.png"))
        plt.close(fig)
    return files


def main(out: Path | None = None, ncount: float = NCOUNT, seed: int = 1) -> list[Path]:
    """Run McStas in every mode with :func:`simulate`, then reduce and draw with :func:`analyze`.

    Parameters
    ----------
    out:
        Folder of the files, as :func:`analyze`.
    ncount, seed:
        As :func:`simulate`.

    Returns
    -------
    :
        The paths of the files.
    """
    simulate(ncount=ncount, seed=seed)
    return analyze(out)


if __name__ == "__main__":
    analyze(out=cli.parse_out(__doc__))  # McStas runs only with main
