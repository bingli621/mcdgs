"""The elastic resolution of T-Rex, from vanadium: the intensity of each sub-pulse against E.

:class:`TRexElasticResolution` is :class:`.TRex` with the vanadium cylinder of :func:`.vanadium`
on ``sample_pos``, with the beamstop and the banana detector. Its control parameters are those
of :class:`.TRex`, so the choppers turn and the beam on the sample is a train of sub-pulses, each
with its own incident energy: unlike the white beam of :mod:`.diagnostics.trex.solid_angles`.
Vanadium scatters incoherently and elastically, in every direction, so every neutron it scatters
has an energy transfer of zero: the width of its elastic line is the energy resolution of the
instrument, and whatever lies away from it comes from the reduction or the instrument, not the
sample. The cylinder scatters only once, into the window of
:func:`.solid_angles.detector_window`, which covers the banana detector.

:func:`main` runs it in two parts, with :class:`.Simulation`, as
:mod:`.diagnostics.trex.elastic_plane` does: from the source to ``sample_pos``, which dumps the
beam there into a beam file, then from that beam through the vanadium to the detector.

The reduction is that of :mod:`.diagnostics.trex.elastic_plane`:

1. :func:`.elastic_plane.subpulse_timing` finds the sub-pulses at the pre-sample and beamstop
   monitors of the first part, with their incident energy and the time they reach the sample.
2. :func:`.read_events` reads the neutrons of the banana detector of the second part, by pixel,
   and :func:`scattered_pixels` leaves out the pixels within ``MIN_TWO_THETA`` of the beam: the
   direct beam that passes the beamstop.
3. :func:`.elastic_plane.normalized_events` keeps the neutrons inside the arrival window of a
   sub-pulse, gives them their energy transfer, and divides each by the solid angle of its pixel
   and by the flux of its sub-pulse.
4. :func:`energy_spectra` sums them into a histogram of the energy transfer of each sub-pulse,
   averaged over the solid angle of the pixels: per incident neutron, per steradian and per meV.

:func:`plot_spectra` draws the spectra, and :func:`elastic_lines` gives the center, the width
and the integral of the elastic line of each sub-pulse.

:func:`main` is :func:`simulate`, which runs McStas, then :func:`analyze`, which reduces and
draws. Run it to do both::

    pixi run python -c "from mcdgs.diagnostics.trex.elastic_resolution_incoherent import main; main()"

Running the module only reduces the runs again, without McStas, with :func:`analyze`::

    pixi run python -m mcdgs.diagnostics.trex.elastic_resolution_incoherent
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipp as sc
from matplotlib.figure import Figure

from mcdgs.core.components import Incoherent
from mcdgs.diagnostics import cli
from mcdgs.diagnostics.trex import elastic_plane as plane
from mcdgs.diagnostics.trex.solid_angles import detector_window
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.mcstas import Simulation, Workspace, read_events
from mcdgs.reduction.binning import bin_edges
from mcdgs.reduction.inelastic import DEFAULT_ENERGY_TRANSFER_RATIO
from mcdgs.samples import vanadium

NAME = "trex_elastic_resolution_incoherent"
TO_SAMPLE = f"{NAME}_run_to_sample"  # the run folder of the first part, source -> sample_pos
FROM_SAMPLE = f"{NAME}_run_from_sample"  # that of the second part, sample_pos -> detector
MIN_TWO_THETA = 5.0  # deg: the pixels closer to the beam see the direct beam past the beamstop
ENERGY_STEP = 0.05  # meV, the bin width of the spectra
ZOOM = 0.1  # the elastic lines span this fraction of the highest Ei on each side of E = 0
LOG_DECADES = plane.LOG_DECADES  # the log scale spans this many decades below the highest bin


# --8<-- [start:instrument]
def sample() -> Incoherent:
    """Return the vanadium cylinder of :func:`.vanadium`, scattering once into the detector.

    It scatters into the window of :func:`.solid_angles.detector_window` only.
    """
    return vanadium(name="sample", order=1, **detector_window())


class TRexElasticResolution(TRex):
    """T-Rex with a vanadium cylinder, from the source to the banana detector.

    The cylinder is that of :func:`sample`.

    Parameters
    ----------
    params:
        Control parameters.
    source_masked:
        As in :class:`.TRex`. True by default: the masked source of
        :func:`.sources.masked_source` spends the rays where the choppers let neutrons through.
    """

    name = "trex_elastic_resolution"  # the McStas instrument name

    def __init__(self, params: Params = DEFAULT_PARAMS, source_masked: bool = True) -> None:
        super().__init__(sample=sample(), params=params, source_masked=source_masked)


def define_instrument() -> TRexElasticResolution:
    """Return the instrument of this diagnostic: :class:`TRexElasticResolution` at its defaults."""
    params = Params()  # the default control parameters
    return TRexElasticResolution(params=params)


# --8<-- [end:instrument]


def two_theta(binned: sc.DataArray) -> sc.Variable:
    """Return the scattering angle of each pixel, from the beam along ``BEAM_DIRECTION``.

    Parameters
    ----------
    binned:
        The neutrons by pixel, as :func:`.read_events`, with the pixel ``position`` and the
        ``sample_position``, in the beam frame of ``sample_pos``.

    Returns
    -------
    :
        The angle between the beam and the direction from the sample to each pixel, in deg,
        along the pixel dimension.
    """
    scattered = binned.coords["position"] - binned.coords["sample_position"]
    beam = sc.vector(plane.BEAM_DIRECTION)
    cosine = sc.dot(scattered, beam) / sc.norm(scattered)
    return sc.acos(sc.values(cosine)).to(unit="deg")


def scattered_pixels(binned: sc.DataArray, min_two_theta: float = MIN_TWO_THETA) -> sc.DataArray:
    """Leave out the pixels close to the beam, which see the direct beam past the beamstop.

    Parameters
    ----------
    binned:
        The neutrons by pixel, as :func:`.read_events`.
    min_two_theta:
        The smallest scattering angle kept, in deg.

    Returns
    -------
    :
        The pixels at ``two_theta`` of ``min_two_theta`` or more, from :func:`two_theta`.
    """
    (dim,) = binned.dims
    keep = two_theta(binned).values >= min_two_theta
    return binned[dim, np.flatnonzero(keep)]


def energy_spectra(
    events: sc.DataArray,
    pulses: sc.DataArray,
    solid_angle: sc.Variable,
    step: float = ENERGY_STEP,
) -> sc.DataArray:
    """Sum the neutrons of each sub-pulse into a histogram of their energy transfer.

    Each neutron of :func:`.elastic_plane.normalized_events` times the solid angle of its pixel
    is the fraction of the incident neutrons of its sub-pulse it stands for. Summed over the
    pixels and divided by their solid angle, it is the intensity averaged over the detector.

    Parameters
    ----------
    events:
        The neutrons of :func:`.elastic_plane.normalized_events`, in 1/sr, with their
        ``subpulse``, ``energy_transfer`` and the ``solid_angle`` of their pixel.
    pulses:
        The sub-pulses of :func:`.elastic_plane.subpulse_timing`.
    solid_angle:
        The solid angle of all the pixels the neutrons come from, in sr.
    step:
        The bin width of the energy transfer, in meV.

    Returns
    -------
    :
        The intensity per incident neutron, per steradian and per meV, along ``subpulse`` and
        ``energy_transfer``, with the bin edges ``energy_transfer`` in meV, symmetric about 0 and
        wide enough for the arrival window of the highest Ei, and the ``incident_energy`` and
        ``incident_wavelength`` of each sub-pulse.
    """
    highest = float(pulses.coords["incident_energy"].to(unit="meV").max().value)
    limit = np.ceil(max(map(abs, DEFAULT_ENERGY_TRANSFER_RATIO)) * highest)
    edges = bin_edges("energy_transfer", -limit, limit, step, unit="meV")
    fractions = events.copy(deep=False)
    fractions.data = events.data * events.coords["solid_angle"]  # of the incident neutrons
    fractions.coords["energy_transfer"] = events.coords["energy_transfer"].to(unit="meV")
    n = pulses.sizes["subpulse"]
    grouped = fractions.group(sc.arange("subpulse", n))
    spectra = grouped.hist(energy_transfer=edges) / solid_angle.to(unit="sr")
    spectra = spectra / sc.scalar(step, unit="meV")
    for name in ("incident_energy", "incident_wavelength"):
        spectra.coords[name] = pulses.coords[name].copy()
    return spectra


def half_height_width(centers: np.ndarray, values: np.ndarray) -> tuple[float, float]:
    """Return where a line peaks and its full width at half its height.

    The half heights are interpolated linearly between the bins on either side of them, on each
    flank of the highest bin.

    Parameters
    ----------
    centers:
        The centers of the bins.
    values:
        The height of each bin.

    Returns
    -------
    :
        The center of the highest bin, and the width, in the unit of ``centers``.

    Examples
    --------
    A triangle 2 wide at its base is 1 wide at half its height:

    >>> centers = np.arange(-1.75, 2.0, 0.5)
    >>> half_height_width(centers, np.array([0.0, 0.0, 0.0, 0.5, 1.0, 0.5, 0.0, 0.0]))
    (0.25, 1.0)
    """
    top = int(np.argmax(values))
    half = values[top] / 2
    below = np.flatnonzero(values[:top] <= half)  # the rising flank
    if below.size:
        i = below[-1]
        rise = float(np.interp(half, values[i : i + 2], centers[i : i + 2]))
    else:
        rise = float(centers[0])
    above = np.flatnonzero(values[top:] <= half)  # the falling flank
    if above.size:
        j = top + above[0]
        fall = float(np.interp(half, values[j - 1 : j + 1][::-1], centers[j - 1 : j + 1][::-1]))
    else:
        fall = float(centers[-1])
    return float(centers[top]), fall - rise


def elastic_lines(spectra: sc.DataArray) -> sc.DataArray:
    """Return the center, the width and the integral of the elastic line of each sub-pulse.

    Parameters
    ----------
    spectra:
        The spectra of :func:`energy_spectra`, along ``subpulse`` and ``energy_transfer``.

    Returns
    -------
    :
        The integral of each spectrum over the energy transfer, per incident neutron and per
        steradian, along ``subpulse``, with the coordinates ``center`` and ``fwhm`` of
        :func:`half_height_width`, in meV, and those of ``spectra`` along ``subpulse``.
    """
    spectra = spectra.transpose(["subpulse", "energy_transfer"])
    edges = spectra.coords["energy_transfer"].to(unit="meV")
    centers = sc.midpoints(edges).values
    found = [half_height_width(centers, values) for values in spectra.values]
    lines = (spectra * (edges[1:] - edges[:-1])).sum("energy_transfer")
    for name, values in zip(("center", "fwhm"), zip(*found)):
        lines.coords[name] = sc.array(dims=["subpulse"], values=list(values), unit="meV")
    return lines


def lines_table(lines: sc.DataArray) -> str:
    """Return the elastic line of each sub-pulse as a Markdown table.

    Parameters
    ----------
    lines:
        The elastic lines of :func:`elastic_lines`.

    Returns
    -------
    :
        One row per sub-pulse: its λi and Ei, the center and the width of its line, the width as
        a fraction of Ei, and the integral of its spectrum.
    """
    rows = []
    for i, (wavelength, ei, center, fwhm, integral) in enumerate(
        zip(
            lines.coords["incident_wavelength"].to(unit="Å").values,
            lines.coords["incident_energy"].to(unit="meV").values,
            lines.coords["center"].values,
            lines.coords["fwhm"].values,
            lines.values,
        )
    ):
        rows.append(
            [
                f"{i}",
                f"{wavelength:.3f} Å",
                f"{ei:.2f} meV",
                f"{center:.3f} meV",
                f"{fwhm:.3f} meV",
                f"{100 * fwhm / ei:.2f} %",
                f"{integral:.3e}",
            ]
        )
    header = ["Sub-pulse", "λi", "Ei", "Center", "FWHM", "FWHM / Ei", "Intensity [1/sr]"]
    return cli.markdown_table(header, rows)


def _subpulse_colors(n: int) -> list[tuple[float, float, float, float]]:
    """Return a color of turbo for each of ``n`` sub-pulses, from the shortest wavelength."""
    cmap = plt.get_cmap("turbo")
    return [cmap(x) for x in np.linspace(0.05, 0.95, n)]


def plot_spectra(
    spectra: sc.DataArray,
    title: str,
    zoom: float = ZOOM,
    window: tuple[float, float] = DEFAULT_ENERGY_TRANSFER_RATIO,
) -> Figure:
    """Draw the spectrum of every sub-pulse against the energy transfer.

    The top panel draws the whole arrival window of every sub-pulse, against the energy transfer
    as a fraction of its Ei, on the log scale of ``LOG_DECADES`` below the highest bin; the
    bottom one each elastic line, against the energy transfer, on a linear scale.

    Parameters
    ----------
    spectra:
        The spectra of :func:`energy_spectra`.
    title:
        The title of the figure.
    zoom:
        The span of the bottom panel on each side of the elastic line, as a fraction of the
        highest Ei.
    window:
        The smallest and the largest ``E / Ei`` of the arrival windows, the span of the top
        panel.

    Returns
    -------
    :
        The figure.
    """
    spectra = spectra.transpose(["subpulse", "energy_transfer"])
    edges = spectra.coords["energy_transfer"].to(unit="meV").values
    centers = (edges[:-1] + edges[1:]) / 2
    energies = spectra.coords["incident_energy"].to(unit="meV").values
    colors = _subpulse_colors(spectra.sizes["subpulse"])
    top = float(spectra.values.max())
    unit = f"Intensity [{spectra.unit}]"

    fig, (full, line) = plt.subplots(2, 1, figsize=(10, 10), layout="constrained")
    fig.suptitle(title)
    for i, (values, ei, color) in enumerate(zip(spectra.values, energies, colors)):
        inside = (centers >= window[0] * ei) & (centers <= window[1] * ei)
        shown = np.ma.masked_less_equal(values[inside], 0)  # empty bins left out
        full.step(centers[inside] / ei, shown, where="mid", color=color, lw=0.8)
        line.plot(centers, values, color=color, lw=1.0, label=f"{i}: Ei = {ei:.1f} meV")
    full.set(
        xlim=window,
        yscale="log",
        ylim=(top / 10**LOG_DECADES, 2 * top),
        xlabel="E / Ei",
        ylabel=unit,
        title=f"The whole window of every sub-pulse, {window[0]:g} Ei to {window[1]:g} Ei".replace(
            "-", "−"
        ),
    )
    reach = zoom * energies.max()
    line.set(
        xlim=(-reach, reach),
        xlabel="Energy transfer E = Ei − Ef [meV]",
        ylabel=unit,
        title="The elastic line of each sub-pulse",
    )
    for ax in (full, line):
        ax.grid(True, alpha=0.3)
    line.legend(fontsize="small", ncols=3)
    return fig


def simulate(ncount: float = 1e9, seed: int = 1) -> tuple[Path, Path]:
    """Run :class:`TRexElasticResolution` in McStas, in two parts split at ``sample_pos``.

    The first part runs from the source to ``sample_pos`` into ``TO_SAMPLE``, with the monitors
    of the incident energy and the beam file; the second from that beam through the vanadium to
    the detector, into ``FROM_SAMPLE``, both in the runs folder of :meth:`.Workspace.find`. Each
    replaces the folder of an earlier run.

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
        The sub-pulses of :func:`.elastic_plane.subpulse_timing`, and the spectra of
        :func:`energy_spectra`, from the pixels of :func:`scattered_pixels`.
    """
    workspace = Workspace.find()
    to_sample = Path(to_sample or workspace.run_folder(TO_SAMPLE))
    from_sample = Path(from_sample or workspace.run_folder(FROM_SAMPLE))
    pulses = plane.subpulse_timing_of_run(to_sample)
    binned = read_events(
        from_sample, plane.DETECTOR, "source", plane.SAMPLE, beam_direction=plane.BEAM_DIRECTION
    )
    binned = scattered_pixels(binned)
    events = plane.normalized_events(binned, pulses)
    solid_angle = plane.pixel_solid_angle(binned).sum()
    return pulses, energy_spectra(events, pulses, solid_angle)


def analyze(
    out: Path | None = None,
    to_sample: str | Path | None = None,
    from_sample: str | Path | None = None,
) -> list[Path]:
    """Reduce the two parts of a McStas run, and save the spectra and their elastic lines.

    It runs no McStas, so it can run again on the same runs; :func:`main` runs them first. The
    files are ``trex_elastic_resolution_incoherent.png``, from :func:`plot_spectra`, and
    ``trex_elastic_resolution_incoherent_lines.txt``, the table of :func:`lines_table`.

    Parameters
    ----------
    out:
        Folder of the files. If ``None``, ``trex_elastic_resolution_incoherent`` in the results
        folder of :meth:`.Workspace.find`.
    to_sample, from_sample:
        The folders of the two parts, as :func:`reduce_runs`.

    Returns
    -------
    :
        The paths of the files.
    """
    pulses, spectra = reduce_runs(to_sample, from_sample)
    folder = cli.output_folder(out, NAME)
    table = lines_table(elastic_lines(spectra))
    print(table, end="")
    files = [cli.save_text(table, folder / f"{NAME}_lines.txt")]
    fig = plot_spectra(
        spectra,
        f"T-Rex, vanadium: the energy transfer of {pulses.sizes['subpulse']} sub-pulses\n"
        f"{define_instrument().running_conditions()}",
    )
    files.append(cli.save_figure(fig, folder / f"{NAME}.png"))
    return files


def main(out: Path | None = None, ncount: float = 1e9, seed: int = 1) -> list[Path]:
    """Run McStas with :func:`simulate`, then reduce and draw it with :func:`analyze`.

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
    to_sample, from_sample = simulate(ncount=ncount, seed=seed)
    return analyze(out, to_sample, from_sample)


if __name__ == "__main__":
    analyze(out=cli.parse_out(__doc__))  # McStas runs only with main, a few minutes
