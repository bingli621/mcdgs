"""A rotation scan of the La2Ni7 crystal on T-Rex, mapping the horizontal plane of reciprocal space.

:class:`TRexWithRotatingLa2Ni7` is :class:`.TRex` with the La2Ni7 crystal of :func:`.la2ni7`,
(2 1 0) in the horizontal plane, on a goniometer that turns it by the instrument parameter
``omega`` about the vertical. Rotating about the vertical sweeps the plane of 2a* + b* and c*
through the horizontal plane of the instrument, so every angle sees more of it.

:func:`simulate` runs only the part from ``sample_pos`` to the detector, once per angle of
``OMEGAS``, from the beam file of :mod:`.diagnostics.trex.elastic_plane`: the beam before the sample does
not depend on the crystal. The crystal repeats each neutron that reaches it ``SCAN_SPLIT``
times, fewer than in :mod:`.diagnostics.trex.elastic_plane`, since every angle lists its own neutrons.

:func:`reduce_angles` reduces each angle as :mod:`.diagnostics.trex.elastic_plane` does, per incident
neutron and per steradian, with the sub-pulses of its run to the sample, turns the momentum
transfer into the frame of the crystal, as it stands at ``omega`` = 0, and keeps the neutrons
near the horizontal plane, ``|Qy| < QY_HALF_WIDTH``, in one (Qx, Qz) map per angle.
:func:`analyze` draws them with :func:`animate`, as a GIF that adds one angle per frame, and
the whole scan, its last frame, with :func:`plot_scan`, both on the log color scale of
:func:`.draw_plane`. :func:`main` is :func:`simulate`, then :func:`analyze`. Run the module
for both::

    pixi run python -m mcdgs.diagnostics.trex.elastic_scan

or reduce the scan again without McStas::

    pixi run python -c "from mcdgs.diagnostics.trex.elastic_scan import analyze; analyze()"
"""

from pathlib import Path

import matplotlib.pyplot as plt
import mcstastox
import numpy as np
import scipp as sc
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.artist import Artist
from matplotlib.figure import Figure

from mcdgs.core.components import Single_crystal
from mcdgs.core.geometry import rotation
from mcdgs.diagnostics.cli import output_folder, parse_out, save_figure
from mcdgs.diagnostics.trex import elastic_plane as plane
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.mcstas.nexus_io import numeric_parameters, read_events
from mcdgs.mcstas.run_folders import default_output_dir
from mcdgs.mcstas.simulation import Simulation
from mcdgs.reduction.binning import bin_edges
from mcdgs.samples import la2ni7

NAME = "trex_elastic_scan"
SCAN = f"{NAME}_omega_from_sample"  # the run folder of the scan, sample_pos -> detector
OMEGAS = np.arange(-30.0, 31.0, 2.0)  # deg, the angles of the scan
IN_PLANE = [2, 1, 0]  # hkl of the reciprocal vector in the horizontal plane at omega = 0
SCAN_SPLIT = 10  # McStas SPLIT on the crystal, at each angle
Q_LIMIT = (
    10.0  # 1/Å: the map spans -Q_LIMIT to Q_LIMIT in Qx and Qz, beyond 2 ki of every sub-pulse
)
Q_STEP = plane.Q_STEP
FRAMES_PER_SECOND = 4  # of the GIF, one frame per angle
QY_HALF_WIDTH = plane.QY_HALF_WIDTH


class TRexWithRotatingLa2Ni7(TRex):
    """T-Rex with a La2Ni7 crystal on a goniometer, turned by ``omega`` about the vertical.

    The crystal is that of :func:`.la2ni7`, with ``IN_PLANE`` in the horizontal plane at
    ``omega`` = 0, split ``SCAN_SPLIT`` times.

    Parameters
    ----------
    params:
        Control parameters.
    source_masked:
        As in :class:`.TRex`. True by default: the masked source of
        :func:`.sources.masked_source` spends the rays where the choppers let neutrons through.
    """

    name = "trex_la2ni7_scan"  # the McStas instrument name

    def __init__(self, params: Params = DEFAULT_PARAMS, source_masked: bool = True) -> None:
        super().__init__(sample=crystal(), params=params, source_masked=source_masked)


def crystal(split: int = SCAN_SPLIT) -> Single_crystal:
    """Return the La2Ni7 crystal of the scan, turned by ``omega`` about the vertical.

    Parameters
    ----------
    split:
        McStas SPLIT on the crystal.

    Returns
    -------
    :
        The crystal, with ``IN_PLANE`` in the horizontal plane at ``omega`` = 0.
    """
    return la2ni7(name="sample", split=split, in_plane=IN_PLANE, rotated=[0, "omega", 0])


def simulate(
    omegas: np.ndarray = OMEGAS, beam_file: str | Path | None = None, seed: int = 1
) -> Path:
    """Run the scan from ``sample_pos`` to the detector, once per angle, from a beam file.

    Parameters
    ----------
    omegas:
        The angles of the crystal, in deg.
    beam_file:
        The beam at ``sample_pos``. If ``None``, that of :mod:`.diagnostics.trex.elastic_plane`, in its
        run to the sample.
    seed:
        Seed of the random numbers; angle ``n`` takes ``seed + n + 1``.

    Returns
    -------
    :
        The scan folder, ``SCAN`` in :func:`.default_output_dir`, with one entry per angle.
    """
    beam_file = beam_file or default_output_dir() / plane.TO_SAMPLE / "beam.mcpl.gz"
    simulation = Simulation(TRexWithRotatingLa2Ni7())  # the default control parameters
    scan_file = simulation.scan_from_sample(
        "omega", omegas, folder_name=SCAN, beam_file=beam_file, seed=seed
    )
    return scan_file.parent


def crystal_frame(events: sc.DataArray, omega: float) -> sc.DataArray:
    """Turn the momentum transfer of neutrons into the frame of a crystal turned by ``omega``.

    The frame is that of the crystal at ``omega`` = 0, the beam frame of ``sample_pos``, so
    every angle of a scan lands its reflections at the same Q.

    Parameters
    ----------
    events:
        The neutrons of :func:`.events_to_q_and_energy`, with ``Qx``, ``Qy`` and ``Qz``.
    omega:
        The angle of the crystal about the vertical, in deg, as McStas ``ROTATED (0, omega, 0)``.

    Returns
    -------
    :
        The neutrons, with ``Qx``, ``Qy`` and ``Qz`` in the frame of the crystal.

    Examples
    --------
    Turned by 90° about the vertical, the crystal's x is the lab's -z:

    >>> events = sc.DataArray(
    ...     sc.ones(dims=["event"], shape=[1]),
    ...     coords={
    ...         name: sc.array(dims=["event"], values=[value], unit="1/Å")
    ...         for name, value in (("Qx", 0.0), ("Qy", 0.0), ("Qz", -1.0))
    ...     },
    ... )
    >>> turned = crystal_frame(events, 90.0)
    >>> [round(float(turned.coords[n].values[0]), 6) + 0.0 for n in ("Qx", "Qy", "Qz")]
    [1.0, 0.0, 0.0]
    """
    names = ("Qx", "Qy", "Qz")
    unit = events.coords["Qx"].unit
    lab = np.stack([events.coords[name].to(unit=unit).values for name in names])
    turned = rotation(0.0, omega, 0.0) @ lab  # McStas: from the lab frame to the crystal's
    events = events.copy(deep=False)
    for name, values in zip(names, turned):
        events.coords[name] = sc.array(dims=events.dims, values=values, unit=unit)
    return events


def reduce_angles(
    scan: str | Path | None = None,
    to_sample: str | Path | None = None,
    step: float = Q_STEP,
    qy_half_width: float = QY_HALF_WIDTH,
) -> sc.DataArray:
    """Reduce every angle of a scan into its own (Qx, Qz) map, in the crystal's frame.

    The sub-pulses come from the run to the sample of :mod:`.diagnostics.trex.elastic_plane`, the same
    for every angle. Each angle is read and reduced on its own, so the scan never needs all its
    neutrons in memory at once.

    Parameters
    ----------
    scan:
        The scan folder of :func:`simulate`. If ``None``, ``SCAN`` in :func:`.default_output_dir`.
    to_sample:
        The run to the sample. If ``None``, that of :mod:`.diagnostics.trex.elastic_plane`.
    step:
        The bin width in Qx and Qz, in 1/Å.
    qy_half_width:
        Half the width of the slice in Qy, in the crystal's frame, in 1/Å.

    Returns
    -------
    :
        The intensity per bin of :func:`.normalized_events`, per incident neutron and per
        steradian, summed over the slice in Qy, along ``omega``, ``Qz`` and ``Qx``, with the
        angle of each map as the coordinate ``omega``, in deg, and the bin edges ``Qz`` and
        ``Qx``.

    Raises
    ------
    ValueError
        If the scan has no angles.
    """
    scan = Path(scan or default_output_dir() / SCAN)
    pulses = plane.subpulse_timing_of_run(to_sample)
    edges = {name: bin_edges(name, -Q_LIMIT, Q_LIMIT, step, unit="1/Å") for name in ("Qz", "Qx")}
    with mcstastox.Read(str(scan)) as data:
        n_entries = data.get_number_of_entries()
    if n_entries == 0:
        raise ValueError(f"{scan}: the scan has no angles")
    maps, omegas = [], []
    for entry in range(1, n_entries + 1):
        with mcstastox.Read(str(scan), entry_number=entry) as data:
            omega = numeric_parameters(data.get_instrument_parameters())["omega"]
        binned = read_events(
            scan,
            plane.DETECTOR,
            "source",
            plane.SAMPLE,
            beam_direction=plane.BEAM_DIRECTION,
            entry=entry,
        )
        events = crystal_frame(plane.normalized_events(binned, pulses), omega)
        qy = events.coords["Qy"].to(unit="1/Å").values
        (dim,) = events.dims
        near = events[dim, np.flatnonzero(np.abs(qy) < qy_half_width)]
        maps.append(near.hist(Qz=edges["Qz"], Qx=edges["Qx"]))
        omegas.append(omega)
        print(f"omega {omega:+6.1f}°: {near.sizes[dim]:9d} neutrons near the plane")
    stack = sc.concat(maps, "omega")
    stack.coords["omega"] = sc.array(dims=["omega"], values=omegas, unit="deg")
    return stack


def crop(plane: sc.DataArray) -> sc.DataArray:
    """Cut a (Qz, Qx) histogram to the bins from the first to the last with neutrons.

    A histogram with dimensions before Qz and Qx, such as ``omega``, is cut to the bins with
    neutrons in any of them.

    Examples
    --------
    >>> plane = sc.DataArray(
    ...     sc.array(dims=["Qz", "Qx"], values=[[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]),
    ...     coords={"Qz": sc.arange("Qz", 4.0), "Qx": sc.arange("Qx", 4.0)},
    ... )
    >>> crop(plane).sizes
    {'Qz': 1, 'Qx': 1}
    """
    values = plane.transpose([*(d for d in plane.dims if d not in ("Qz", "Qx")), "Qz", "Qx"])
    filled = (values.values > 0).reshape(-1, *values.values.shape[-2:]).any(axis=0)
    rows = np.flatnonzero(filled.any(axis=1))
    columns = np.flatnonzero(filled.any(axis=0))
    if rows.size == 0:
        return plane
    return plane["Qz", rows[0] : rows[-1] + 1]["Qx", columns[0] : columns[-1] + 1]


def angles_summed(omegas: np.ndarray, n: int) -> str:
    """Describe the first ``n`` angles of a scan, such as in the title of their sum.

    Examples
    --------
    >>> angles_summed(np.array([-30.0, -28.0, -26.0]), 2)
    'omega from -30.0° to -28.0°, 2 angles summed'
    """
    return f"omega from {omegas[0]:+.1f}° to {omegas[n - 1]:+.1f}°, {n} angles summed"


def plot_scan(maps: sc.DataArray, title: str) -> Figure:
    """Draw the whole scan: the sum of the maps of every angle, as :func:`.draw_plane`.

    Parameters
    ----------
    maps:
        The maps of :func:`reduce_angles`, along ``omega``, in the order of the scan.
    title:
        The title of the figure; the angles summed are added below it.

    Returns
    -------
    :
        The figure, Qx along ``IN_PLANE`` and Qz along (0 0 1).
    """
    hkl = " ".join(str(i) for i in IN_PLANE)
    omegas = maps.coords["omega"].to(unit="deg").values
    return plane.draw_plane(
        maps.sum("omega"),
        f"{title}\n{angles_summed(omegas, len(omegas))}",
        xlabel=f"Qx [1/Å], along ({hkl})",
        ylabel="Qz [1/Å], along (0 0 1)",
    )


def animate(maps: sc.DataArray, path: Path, title: str, fps: int = FRAMES_PER_SECOND) -> Path:
    """Save the scan as a GIF that adds one angle per frame, with the angles so far in its title.

    Frame ``n`` is the sum of the maps of the first ``n + 1`` angles, so the map fills in as the
    crystal turns, and the last frame is the whole scan, :func:`plot_scan`. Every frame has
    its color scale, so that the frames compare.

    Parameters
    ----------
    maps:
        The maps of :func:`reduce_angles`, along ``omega``, in the order of the scan.
    path:
        The GIF to write.
    title:
        The title of every frame; the angles so far are added below it.
    fps:
        Frames per second.

    Returns
    -------
    :
        ``path``.
    """
    omegas = maps.coords["omega"].to(unit="deg").values
    values = maps.transpose(["omega", "Qz", "Qx"]).values
    summed = np.cumsum(values, axis=0)  # frame n: the angles up to n
    fig = plot_scan(maps, title)  # the axes and the color scale of the whole scan
    (mesh,) = fig.axes[0].collections

    def frame(n: int) -> list[Artist]:
        mesh.set_array(np.ma.masked_less_equal(summed[n], 0))
        fig.suptitle(f"{title}\n{angles_summed(omegas, n + 1)}")
        return [mesh]

    animation = FuncAnimation(fig, frame, frames=len(omegas))
    animation.save(path, writer=PillowWriter(fps=fps), dpi=100)
    plt.close(fig)
    return path


def analyze(
    out: Path | None = None,
    scan: str | Path | None = None,
    to_sample: str | Path | None = None,
) -> list[Path]:
    """Reduce a scan and save it as ``trex_elastic_scan.gif`` and ``trex_elastic_scan.png``.

    The GIF, from :func:`animate`, adds one angle per frame; the figure, from
    :func:`plot_scan`, is its last frame, the whole scan. It runs no McStas, so it can run
    again on the same scan.

    Parameters
    ----------
    out:
        Folder of the files. If ``None``, the scan folder.
    scan, to_sample:
        As :func:`reduce_angles`.

    Returns
    -------
    :
        The paths of the GIF and of the figure.
    """
    scan = Path(scan or default_output_dir() / SCAN)
    maps = crop(reduce_angles(scan, to_sample))
    hkl = " ".join(str(i) for i in IN_PLANE)
    title = (
        f"T-Rex, La2Ni7 ({hkl}) in plane\n{TRexWithRotatingLa2Ni7().running_conditions()}\n"
        f"|Qy| < {QY_HALF_WIDTH:g} 1/Å, in the frame of the crystal at omega = 0"
    )
    folder = output_folder(out, scan)
    gif = animate(maps, folder / f"{NAME}.gif", title)
    print(f"Saved {gif}")
    return [gif, save_figure(plot_scan(maps, title), folder / f"{NAME}.png")]


def main(out: Path | None = None) -> list[Path]:
    """Run the scan with :func:`simulate`, then reduce and draw it with :func:`analyze`.

    The beam file must exist: run :mod:`.diagnostics.trex.elastic_plane` first.

    Parameters
    ----------
    out:
        Folder of the files, as :func:`analyze`.

    Returns
    -------
    :
        The paths of the GIF and of the figure.
    """
    return analyze(out, simulate())


if __name__ == "__main__":
    main(out=parse_out(__doc__))
