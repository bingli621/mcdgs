"""The S(Q, E) of a powder magnon in every sub-pulse of CSPEC: the resolution across the train.

:class:`CSPECMagnonPowder` is :class:`.CSPEC` with the powder of :func:`.magnon_powder` on
``sample_pos``, and the detector panels, as
:mod:`.diagnostics.cspec.elastic_resolution_incoherent_short` has its vanadium: the same
cylinder, the same choppers, so the same 14 sub-pulses, each with its own incident energy. It follows :mod:`.diagnostics.trex.powder_magnon`, with a magnon in reach of
the sub-pulses of CSPEC, which have about a third of the Ei of those of T-Rex: the same cell of
:data:`.MAGNON_A`, so the same zone centers in |Q|, with a gap of ``MAGNON_GAP`` and a top of
``MAGNON_TOP``, and a line of ``MAGNON_WIDTH``, far narrower than the resolution of CSPEC. The
powder scatters once, within ``d_phi`` of the horizontal plane: the height of the window of
:func:`.elastic_resolution_incoherent_short.detector_window`, which covers the panels.

:func:`simulate` writes the S(|Q|, E) file of the powder into the builds folder of the workspace,
:func:`sqw_path`, then runs in two parts, with :class:`.Simulation`, as the vanadium does: from
the source to ``sample_pos``, then from the beam there through the powder to the detector, into
the folders of :func:`run_names`.

The reduction is that of the vanadium
(:func:`.elastic_resolution_incoherent_short.reduce_runs`): the sub-pulses at ``mc_monitor``, every pixel of the panels, the neutrons of each sub-pulse normalized
by the solid angle of their pixel and by its flux, and :func:`.q_e_maps`, S(Q, E) normalized by
the coverage of each bin. :func:`.plot_magnon_maps` draws the 14 maps in one figure, over
their whole arrival windows, energy gain and energy loss, on a color scale up to the
``MAP_PERCENTILE`` percentile of the bins, with the edges of the band of the model,
:func:`.magnon_band`, over them.

:func:`main` is :func:`simulate`, which runs McStas, then :func:`analyze`, which reduces and
draws. Run it to do both::

    pixi run python -c "from mcdgs.diagnostics.cspec.powder_magnon import main; main()"

Running the module only reduces the runs again, without McStas, with :func:`analyze`::

    pixi run python -m mcdgs.diagnostics.cspec.powder_magnon
"""

from pathlib import Path

import matplotlib.pyplot as plt
import scipp as sc

from mcdgs.core.components import Isotropic_Sqw
from mcdgs.diagnostics import cli
from mcdgs.diagnostics import elastic_resolution as resolution
from mcdgs.diagnostics.cspec import elastic_resolution_incoherent_short as vanadium
from mcdgs.diagnostics.cspec import sample_runs as runs
from mcdgs.diagnostics.powder_magnon import plot_magnon_maps
from mcdgs.instruments.cspec.model import CSPEC, DEFAULT_PARAMS, Params
from mcdgs.mcstas import Simulation, Workspace
from mcdgs.reduction.inelastic import normalized_events
from mcdgs.samples import magnon_powder, write_magnon_sqw

NAME = "cspec_magnon_powder"
SQW_NAME = "cspec_magnon_powder.sqw"  # in the builds folder of the workspace; not that of T-Rex
# the magnon: a gap of MAGNON_GAP at Q = 0, up to MAGNON_TOP at the corner of the zone, within
# the 0.8 Ei that the arrival windows of the sub-pulses reach, 12 to 4.6 meV
MAGNON_GAP = 1.0  # meV
MAGNON_TOP = 8.0  # meV
MAGNON_JS = (MAGNON_TOP - MAGNON_GAP) / 24  # meV: the band is 24 JS wide
MAGNON_WIDTH = 0.01  # meV, the σ of the line: far below the narrowest FWHM of CSPEC, 0.1 meV
SQW_E_STEP = 0.005  # meV, the spacing of E of the table, finer than the line
SQW_Q_MAX = 6.0  # 1/Å, past the largest |Q| of CSPEC, 2 ki at the highest Ei
P_INTERACT = 0.9  # the fraction of the neutrons made to scatter in the powder
SAMPLE_SPLIT = 4  # each neutron of the beam goes through the powder this many times
NCOUNT = 5e9  # McStas rays of the first part, from the masked source
ENERGY_STEP = 0.1  # meV, the bin width in E of the maps
Q_STEP = 0.05  # 1/Å, the bin width in |Q| of the maps
# the top of the color scale, as a percentile of the bins with neutrons: lower than for vanadium,
# so that the band shows in the upper colors, not only the hottest bins of the last sub-pulse
MAP_PERCENTILE = 98.0


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


def write_sqw(path: str | Path) -> Path:
    """Write the S(|Q|, E) of the magnon of this diagnostic with :func:`.write_magnon_sqw`.

    The magnon has ``MAGNON_JS``, ``MAGNON_GAP`` and ``MAGNON_WIDTH``; the table runs to
    ``SQW_Q_MAX`` in |Q|, in steps of ``SQW_E_STEP`` in E.

    Parameters
    ----------
    path:
        The file to write.

    Returns
    -------
    :
        The path of the file.
    """
    return write_magnon_sqw(
        path,
        js=MAGNON_JS,
        gap=MAGNON_GAP,
        width=MAGNON_WIDTH,
        q_max=SQW_Q_MAX,
        e_step=SQW_E_STEP,
    )


# --8<-- [start:instrument]
def sample(sqw_file: str | Path | None = None) -> Isotropic_Sqw:
    """Return the powder of :func:`.magnon_powder`, scattering once towards the detector.

    It scatters within ``d_phi`` of the horizontal plane, the height of the window of
    :func:`.elastic_resolution_incoherent_short.detector_window`, and makes ``P_INTERACT`` of the
    neutrons scatter. Each neutron goes through it ``SAMPLE_SPLIT`` times (McStas SPLIT), to
    fill S(Q, E) more smoothly.

    Parameters
    ----------
    sqw_file:
        Its S(|Q|, E) file, from :func:`write_sqw`. If ``None``, that of :func:`sqw_path`.
    """
    return magnon_powder(
        sqw_path() if sqw_file is None else sqw_file,
        name="sample",
        order=1,
        d_phi=vanadium.detector_window()["focus_ah"],
        p_interact=P_INTERACT,
        split=SAMPLE_SPLIT,
    )


class CSPECMagnonPowder(CSPEC):
    """CSPEC with a powder of a magnon, from the source to the detector panels.

    The powder is that of :func:`sample`. ``mc_monitor`` is :func:`.sample_runs.timing_monitor`,
    with finer time bins, to time each sub-pulse.

    Parameters
    ----------
    params:
        Control parameters.
    source_masked:
        As in :class:`.CSPEC`. True by default: the masked source of
        :func:`.sources.masked_source` spends the rays where the choppers let neutrons through.
    sqw_file:
        As :func:`sample`.
    """

    name = "cspec_magnon_powder"  # the McStas instrument name

    def __init__(
        self,
        params: Params = DEFAULT_PARAMS,
        source_masked: bool = True,
        sqw_file: str | Path | None = None,
    ) -> None:
        super().__init__(sample=sample(sqw_file), params=params, source_masked=source_masked)
        self.mc_monitor = runs.timing_monitor(params)  # in its place, with finer time bins


def define_instrument() -> CSPECMagnonPowder:
    """Return the instrument of this diagnostic: :class:`CSPECMagnonPowder` at its defaults."""
    params = Params()  # the default control parameters
    return CSPECMagnonPowder(params=params)


# --8<-- [end:instrument]


def run_names() -> tuple[str, str]:
    """Return the run folders: the part to ``sample_pos``, and the part from it.

    Examples
    --------
    >>> run_names()[1]
    'cspec_magnon_powder_run_from_sample'
    """
    return f"{NAME}_run_to_sample", f"{NAME}_run_from_sample"


def simulate(ncount: float = NCOUNT, seed: int = 1) -> tuple[Path, Path]:
    """Write the S(|Q|, E) of the powder, then run :class:`CSPECMagnonPowder` in McStas.

    The file is that of :func:`write_sqw`, at :func:`sqw_path`: McStas reads it when the
    instrument starts, in both parts. The first part runs from the source to ``sample_pos``, the
    second from that beam through the powder to the detector panels, both into the folders of
    :func:`run_names`. Each replaces the folder of an earlier run.

    Parameters
    ----------
    ncount:
        Number of McStas rays of the first part, from the masked source. The second part runs
        every neutron of the beam file.
    seed:
        Seed of the random numbers of the first part; the second part takes ``seed + 1``.

    Returns
    -------
    :
        The folders of the first and the second part.
    """
    write_sqw(sqw_path())
    to_name, from_name = run_names()
    sim = Simulation(define_instrument())
    to_sample = sim.run_to_sample(folder_name=to_name, ncount=ncount, seed=seed)
    # sample_pos -> detector, from the beam file, in the frame of sample_pos: the beam frame
    from_sample = sim.run_from_sample(folder_name=from_name, seed=seed + 1)
    return to_sample.folder, from_sample.folder


def reduce_runs(
    to_sample: str | Path | None = None, from_sample: str | Path | None = None
) -> tuple[sc.DataArray, sc.DataArray]:
    """Reduce the two parts of a McStas run of :func:`simulate` to S(Q, E), without McStas.

    As :func:`.elastic_resolution_incoherent_short.reduce_runs`, with the bins ``ENERGY_STEP`` and
    ``Q_STEP``. Every pixel is kept: the direct beam passes between the two panels closest to it.

    Parameters
    ----------
    to_sample, from_sample:
        The folders of the first and the second part. If ``None``, those of :func:`run_names`
        in the runs folder of :meth:`.Workspace.find`, where :func:`simulate` writes.

    Returns
    -------
    :
        The sub-pulses of :func:`.sample_runs.subpulse_timing`, and the maps of
        :func:`.q_e_maps`.
    """
    workspace = Workspace.find()
    to_name, from_name = run_names()
    to_sample = Path(to_sample or workspace.run_folder(to_name))
    from_sample = Path(from_sample or workspace.run_folder(from_name))
    pulses = runs.subpulse_timing_of_run(to_sample, define_instrument())  # as simulate runs it
    binned = runs.read_detector(from_sample)
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


def analyze(
    out: Path | None = None,
    to_sample: str | Path | None = None,
    from_sample: str | Path | None = None,
) -> Path:
    """Reduce the runs and draw the S(Q, E) of all the sub-pulses.

    It runs no McStas, so it can run again on the same runs; :func:`main` runs them first. The
    figure is ``cspec_magnon_powder_sqe.png``, from :func:`.plot_magnon_maps`.

    Parameters
    ----------
    out:
        Folder of the figure. If ``None``, ``cspec_magnon_powder`` in the results folder of
        :meth:`.Workspace.find`.
    to_sample, from_sample:
        The folders of the two parts, as :func:`reduce_runs`.

    Returns
    -------
    :
        The path of the figure.
    """
    pulses, maps = reduce_runs(to_sample, from_sample)
    fig = plot_magnon_maps(
        maps,
        f"CSPEC, a powder magnon: S(Q, E) of each of the {pulses.sizes['subpulse']} "
        f"sub-pulses\n{define_instrument().running_conditions()}",
        js=MAGNON_JS,
        gap=MAGNON_GAP,
        percentile=MAP_PERCENTILE,
    )
    path = cli.save_figure(fig, cli.output_folder(out, NAME) / f"{NAME}_sqe.png")
    plt.close(fig)
    return path


def main(out: Path | None = None, ncount: float = NCOUNT, seed: int = 1) -> Path:
    """Run McStas with :func:`simulate`, then reduce and draw with :func:`analyze`.

    Parameters
    ----------
    out:
        Folder of the figure, as :func:`analyze`.
    ncount, seed:
        As :func:`simulate`.

    Returns
    -------
    :
        The path of the figure.
    """
    to_sample, from_sample = simulate(ncount=ncount, seed=seed)
    return analyze(out, to_sample, from_sample)


if __name__ == "__main__":
    analyze(out=cli.parse_out(__doc__))  # McStas runs only with main
