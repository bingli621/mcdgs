"""Diffraction from a La2Ni7 crystal on T-Rex: the horizontal scattering plane, in (Qx, Qz).

Every neutron in the arrival window of its sub-pulse counts, whatever its energy transfer: the
inelastic scattering is weak next to the Bragg peaks, so the plane is a diffraction pattern.

:class:`TRexWithLa2Ni7` is :class:`.TRex` with the La2Ni7 crystal of :func:`.la2ni7` on
``sample_pos``, with the beamstop and the banana detector, which lists every neutron reaching
it. Its control parameters are those of :class:`.TRex`. Few neutrons reach a 5 mm crystal, so
McStas repeats each one that does ``SAMPLE_SPLIT`` times, with its weight divided among them.

:func:`main` runs it in two parts, with :class:`.Simulation`: from the source to ``sample_pos``,
which dumps the beam there into a beam file, then from that beam through the crystal to the
detector. The second part can run again from the same beam, without the long first part.

The reduction is that of :mod:`.diagnostics.trex.sample_runs`, with :mod:`mcdgs.reduction`:

1. :func:`.sample_runs.subpulse_timing` finds the sub-pulses at the pre-sample and beamstop
   monitors of the first part, as :mod:`.diagnostics.trex.incident_energy` does, with their
   incident energy and the time they reach the sample.
2. :func:`.sample_runs.read_detector` reads the neutrons of the banana detector of the second
   part, by pixel, in the frame of ``sample_pos``, which is the beam frame: the origin is the
   sample, z the axis of ``arm_mc``, straight along the beam from the pre-sample monitor to the
   sample, y up and x to the left of the beam.
3. :func:`.events_to_q_and_energy`, within :func:`.normalized_events`, gives each pixel and
   sub-pulse the arrival window of :func:`.arrival_window`, keeps the neutrons that
   :func:`.assign_subpulse` finds inside the window of a sub-pulse, whatever their energy
   transfer within it, and converts them with ``QE_GRAPH`` to the momentum transfer (Qx, Qy,
   Qz).
4. :func:`.normalized_events` divides each neutron by the solid angle of its pixel, from
   :func:`.sample_runs.pixel_solid_angle`, and by the incoming flux of its sub-pulse, its
   intensity at the pre-sample monitor: per incident neutron and per steradian.

:func:`plot_plane` draws the intensity in (Qx, Qz) of every sub-pulse, summed over a slice of
Qy near the horizontal plane, ``|Qy| < QY_HALF_WIDTH``, on a log color scale. Each sub-pulse has
its own incident wavevector ki, so its elastic neutrons lie near its own circle, of radius ki
through the origin, where the detector covers it.

:func:`main` is :func:`simulate`, which runs McStas, then :func:`analyze`, which reduces and
draws. Run it to do both::

    pixi run python -c "from mcdgs.diagnostics.trex.diffraction import main; main()"

Running the module only reduces the runs again, without McStas, with :func:`analyze`;
:func:`reduce_runs` gives the sub-pulses and the neutrons alone::

    pixi run python -m mcdgs.diagnostics.trex.diffraction
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipp as sc
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure

from mcdgs.diagnostics import cli
from mcdgs.diagnostics.trex import sample_runs as runs
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.mcstas import Simulation, Workspace
from mcdgs.reduction.binning import bin_edges, bins_title
from mcdgs.reduction.inelastic import DEFAULT_ENERGY_TRANSFER_RATIO, normalized_events
from mcdgs.samples import la2ni7

NAME = "trex_diffraction"
TO_SAMPLE = f"{NAME}_run_to_sample"  # the run folder of the first part, source -> sample_pos
FROM_SAMPLE = f"{NAME}_run_from_sample"  # that of the second part, sample_pos -> detector
Q_STEP = 0.02  # 1/Å, the bin width in Qx and Qz
LOG_DECADES = 4  # the color scale spans this many decades below the highest bin
QY_HALF_WIDTH = 0.5  # 1/Å: the plane keeps the neutrons with |Qy| below this, near horizontal
DIRECT_BEAM_RADIUS = 0.2  # 1/Å: the direct beam, far above the rest, is left out of the scale
SAMPLE_SPLIT = 10  # McStas SPLIT on the crystal: each neutron reaching it is repeated this often


# --8<-- [start:instrument]
class TRexWithLa2Ni7(TRex):
    """T-Rex with a La2Ni7 crystal, from the source to the banana detector.

    The crystal is that of :func:`.la2ni7`, split ``SAMPLE_SPLIT`` times.

    Parameters
    ----------
    params:
        Control parameters.
    source_masked:
        As in :class:`.TRex`. True by default: the masked source of
        :func:`.sources.masked_source` spends the rays where the choppers let neutrons through.
    """

    name = "trex_la2ni7"  # the McStas instrument name

    def __init__(self, params: Params = DEFAULT_PARAMS, source_masked: bool = True) -> None:
        sample = la2ni7(name="sample", split=SAMPLE_SPLIT)
        super().__init__(sample=sample, params=params, source_masked=source_masked)


def define_instrument() -> TRexWithLa2Ni7:
    """Return the instrument of this diagnostic: :class:`TRexWithLa2Ni7` at its defaults."""
    params = Params()  # the default control parameters
    return TRexWithLa2Ni7(params=params)


# --8<-- [end:instrument]


def print_subpulse_table(pulses: sc.DataArray, events: sc.DataArray) -> None:
    """Print the incident energy, ki and time at the sample of each sub-pulse.

    Parameters
    ----------
    pulses:
        The sub-pulses of :func:`.sample_runs.subpulse_timing`.
    events:
        The neutrons of :func:`.normalized_events`.
    """
    print(
        "  sub-pulse  t_sample [ms]  λi [Å]  Ei [meV]  ki [1/Å]  flux [n/s]  neutrons"
        "  on the detector"
    )
    subpulse = events.coords["subpulse"].values
    n = pulses.sizes["subpulse"]
    counts = np.bincount(subpulse, minlength=n)
    # the fraction of the incident neutrons that reach the detector, in a window
    fractions = np.bincount(
        subpulse, events.values * events.coords["solid_angle"].values, minlength=n
    )
    for i, (t, wavelength, ei, flux, count, fraction) in enumerate(
        zip(
            pulses.coords["time_at_sample"].to(unit="ms").values,
            pulses.coords["incident_wavelength"].to(unit="Å").values,
            pulses.coords["incident_energy"].to(unit="meV").values,
            pulses.values,
            counts,
            fractions,
        )
    ):
        k = 2 * np.pi / wavelength
        print(
            f"  {i:9d}  {t:13.3f}  {wavelength:6.4f}  {ei:8.3f}  {k:8.4f}  {flux:10.4g}"
            f"  {count:8d}  {fraction:15.4g}"
        )


def plot_plane(
    events: sc.DataArray,
    title: str,
    step: float = Q_STEP,
    qy_half_width: float = QY_HALF_WIDTH,
) -> Figure:
    """Draw the intensity in (Qx, Qz) of every sub-pulse, near the horizontal plane.

    Only the neutrons with ``|Qy| < qy_half_width`` count, summed over Qy: the detector also
    sees reflections well above and below the plane, which a sum over all Qy would project
    onto it.

    It is drawn on the log color scale of :func:`draw_plane`, with the bins of
    :func:`plot_bins` under the title.

    Parameters
    ----------
    events:
        The neutrons of :func:`.events_to_q_and_energy`.
    title:
        The title of the figure.
    step:
        The bin width in Qx and Qz, in 1/Å.
    qy_half_width:
        Half the width of the slice in Qy, in 1/Å.

    Returns
    -------
    :
        The figure.
    """
    qy = events.coords["Qy"].to(unit="1/Å").values
    (dim,) = events.dims
    events = events[dim, np.flatnonzero(np.abs(qy) < qy_half_width)]
    edges = {}
    for name in ("Qx", "Qz"):
        values = events.coords[name].to(unit="1/Å").values
        start, stop = np.floor(values.min() * 10) / 10, np.ceil(values.max() * 10) / 10
        edges[name] = bin_edges(name, start, stop, step, unit="1/Å")
    plane = events.hist(Qz=edges["Qz"], Qx=edges["Qx"])  # summed over the slice in Qy
    return draw_plane(plane, f"{title}\n{plot_bins(plane, qy_half_width)}")


def plot_bins(
    plane: sc.DataArray,
    qy_half_width: float = QY_HALF_WIDTH,
    window: tuple[float, float] = DEFAULT_ENERGY_TRANSFER_RATIO,
) -> str:
    """Write the bins of a (Qz, Qx) plane for its title, with :func:`.bins_title`.

    Qx and Qz are binned; the plane is summed over the slice in Qy and over the arrival window of
    each sub-pulse in E / Ei.

    Parameters
    ----------
    plane:
        The histogram, with the bin edges ``Qx`` and ``Qz``.
    qy_half_width:
        Half the width of the slice in Qy, in 1/Å.
    window:
        The smallest and the largest ``E / Ei`` of the arrival windows, as
        :func:`.arrival_window`.

    Returns
    -------
    :
        The bins, as ``Qx=(start, stop, step) [1/Å]``.
    """
    return bins_title(
        {
            "Qx": plane.coords["Qx"],
            "Qz": plane.coords["Qz"],
            "Qy": bin_edges("Qy", -qy_half_width, qy_half_width, unit="1/Å"),
            "E/Ei": bin_edges("E/Ei", *window),
        }
    )


def beyond_direct_beam_max(
    plane: sc.DataArray, direct_beam_radius: float = DIRECT_BEAM_RADIUS
) -> float:
    """Return the highest bin of a (Qz, Qx) histogram beyond the direct beam.

    Parameters
    ----------
    plane:
        The histogram, with the bin edges ``Qx`` and ``Qz``, in 1/Å, and any dimensions before
        them.
    direct_beam_radius:
        The |Q| of the direct beam, in 1/Å, left out.

    Returns
    -------
    :
        The highest value of the bins with ``|Q| > direct_beam_radius``.
    """
    qx, qz = (sc.midpoints(plane.coords[name]).to(unit="1/Å").values for name in ("Qx", "Qz"))
    outside = np.hypot(*np.meshgrid(qz, qx, indexing="ij")) > direct_beam_radius
    return float(plane.values[..., outside].max())


def draw_plane(
    plane: sc.DataArray,
    title: str,
    xlabel: str = "Qx [1/Å]",
    ylabel: str = "Qz [1/Å], along the beam",
    direct_beam_radius: float = DIRECT_BEAM_RADIUS,
    vmax: float | None = None,
) -> Figure:
    """Draw a histogram in (Qx, Qz) on a log color scale.

    The scale is turbo, ``LOG_DECADES`` deep below the highest bin beyond the direct beam,
    ``|Q| > direct_beam_radius``, with the empty bins in white: the direct beam that passes the
    beamstop, a few bins far above all others, would leave the rest at the bottom of the scale,
    so it takes the color of the top instead. The axes have the same scale, so that circles look
    round.

    Parameters
    ----------
    plane:
        The histogram, with the bin edges ``Qx`` and ``Qz``, in 1/Å, evenly spaced.
    title:
        The title of the figure.
    xlabel, ylabel:
        The labels of the axes.
    direct_beam_radius:
        The |Q| of the direct beam, in 1/Å, that the color scale leaves out.
    vmax:
        The top of the color scale. If ``None``, the highest bin beyond the direct beam.

    Returns
    -------
    :
        The figure.
    """
    step = float(np.diff(plane.coords["Qx"].values)[0])
    if vmax is None:
        vmax = beyond_direct_beam_max(plane, direct_beam_radius)

    fig, ax = plt.subplots(figsize=(9, 8), layout="constrained")
    fig.suptitle(title)
    mesh = ax.pcolormesh(
        plane.coords["Qx"].values,
        plane.coords["Qz"].values,
        np.ma.masked_less_equal(plane.values, 0),  # empty bins in white
        cmap=plt.get_cmap("turbo").with_extremes(bad="white"),
        norm=LogNorm(vmin=vmax / 10**LOG_DECADES, vmax=vmax, clip=True),
    )
    ax.set(xlabel=xlabel, ylabel=ylabel, aspect="equal")
    ax.grid(True, alpha=0.3)
    fig.colorbar(
        mesh,
        ax=ax,
        shrink=0.8,
        label=f"Intensity [{plane.unit} per ({step:g} 1/Å)²]",
        extend="max",
    )
    return fig


def simulate(ncount: float = 1e9, seed: int = 1) -> tuple[Path, Path]:
    """Run :class:`TRexWithLa2Ni7` in McStas, in two parts split at ``sample_pos``.

    The first part runs from the source to ``sample_pos`` into ``TO_SAMPLE``, with the monitors
    of the incident energy and the beam file; the second from that beam through the crystal to
    the detector, into ``FROM_SAMPLE``, both in the runs folder of :meth:`.Workspace.find`. Each
    replaces the folder of an earlier run.

    Parameters
    ----------
    ncount:
        Number of McStas rays of the first part. With the masked source nearly every ray reaches
        the sample: 1e8 rays leave about as many neutrons in the beam file as 1e10 of the plain
        source, and 1e9 about ten times as many. The second part, and every angle of
        :mod:`.diagnostics.trex.diffraction_scan`, runs every neutron of the beam file.
    seed:
        Seed of the random numbers of the first part; the second part takes ``seed + 1``.

    Returns
    -------
    :
        The folders of the first and the second part.
    """
    sim = Simulation(define_instrument())
    to_sample = sim.run_to_sample(folder_name=TO_SAMPLE, ncount=ncount, seed=seed)
    # sample_pos -> detector, from the beam file, in the frame of sample_pos: the beam frame
    from_sample = sim.run_from_sample(folder_name=FROM_SAMPLE, seed=seed + 1)
    return to_sample.folder, from_sample.folder


def reduce_runs(
    to_sample: str | Path | None = None, from_sample: str | Path | None = None
) -> tuple[sc.DataArray, sc.DataArray]:
    """Reduce the two parts of a McStas run of :func:`simulate`, without running McStas.

    Parameters
    ----------
    to_sample, from_sample:
        The folders of the first and the second part. If ``None``, ``TO_SAMPLE`` and
        ``FROM_SAMPLE`` in the runs folder of :meth:`.Workspace.find`, where :func:`simulate`
        writes.

    Returns
    -------
    :
        The sub-pulses of :func:`.sample_runs.subpulse_timing`, and the neutrons of
        :func:`.normalized_events`.
    """
    workspace = Workspace.find()
    to_sample = Path(to_sample or workspace.run_folder(TO_SAMPLE))
    from_sample = Path(from_sample or workspace.run_folder(FROM_SAMPLE))
    # the first part: the incident energy and the time at the sample of each sub-pulse
    pulses = runs.subpulse_timing_of_run(to_sample, define_instrument())  # as simulate runs it
    binned = runs.read_detector(from_sample)
    return pulses, normalized_events(binned, pulses, runs.pixel_solid_angle(binned))


def analyze(
    out: Path | None = None,
    to_sample: str | Path | None = None,
    from_sample: str | Path | None = None,
) -> Path:
    """Reduce the two parts of a McStas run, print the sub-pulses and save the (Qx, Qz) plane.

    It runs no McStas, so it can run again on the same runs; :func:`main` runs them first. The
    figure is ``trex_diffraction.png``, from :func:`plot_plane`.

    Parameters
    ----------
    out:
        Folder of the figure. If ``None``, ``trex_diffraction`` in the results folder of
        :meth:`.Workspace.find`.
    to_sample, from_sample:
        The folders of the two parts, as :func:`reduce_runs`.

    Returns
    -------
    :
        The path of the figure.
    """
    from_sample = Path(from_sample or Workspace.find().run_folder(FROM_SAMPLE))
    pulses, events = reduce_runs(to_sample, from_sample)
    print_subpulse_table(pulses, events)

    fig = plot_plane(
        events,
        f"T-Rex, La2Ni7: (Qx, Qz) of {pulses.sizes['subpulse']} sub-pulses\n"
        f"{define_instrument().running_conditions()}",
    )
    return cli.save_figure(fig, cli.output_folder(out, NAME) / f"{NAME}.png")


def main(out: Path | None = None, ncount: float = 1e9, seed: int = 1) -> Path:
    """Run McStas with :func:`simulate`, then reduce and draw it with :func:`analyze`.

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
    analyze(out=cli.parse_out(__doc__))  # McStas runs only with main, about two minutes
