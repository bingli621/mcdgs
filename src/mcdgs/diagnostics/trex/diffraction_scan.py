"""A rotation scan of the La2Ni7 crystal on T-Rex, mapping the horizontal plane of reciprocal space.

:class:`TRexWithRotatingLa2Ni7` is :class:`.TRex` with the La2Ni7 crystal of :func:`.la2ni7`,
(2 1 0) in the horizontal plane, on a goniometer that turns it by the instrument parameter
``omega`` about the vertical. Rotating about the vertical sweeps the plane of 2a* + b* and c*
through the horizontal plane of the instrument, so every angle sees more of it.

:func:`simulate` runs only the part from ``sample_pos`` to the detector, once per angle of
``OMEGAS``, from the beam file of :mod:`.diagnostics.trex.diffraction`: the beam before the sample does
not depend on the crystal. The crystal repeats each neutron that reaches it ``SCAN_SPLIT``
times, fewer than in :mod:`.diagnostics.trex.diffraction`, since every angle lists its own neutrons.

:func:`reduce_angles` reduces each angle with :func:`.diffraction_scan.reduce_scan`, as
:mod:`.diagnostics.trex.diffraction` does, per incident neutron and per steradian, with the
sub-pulses of its run to the sample, turns the momentum transfer into the frame of the crystal,
as it stands at ``omega`` = 0, and keeps the neutrons near the horizontal plane,
``|Qy| < QY_HALF_WIDTH``, in one (Qx, Qz) map per angle. :func:`analyze` draws them with
:func:`.diffraction_scan.animate`, as a GIF that adds one angle per frame, and the whole scan,
its last frame, with :func:`.diffraction_scan.plot_scan`, both on the log color scale of
:func:`.draw_plane`. :func:`main` is :func:`simulate`, then :func:`analyze`. Run the module
for both::

    pixi run python -m mcdgs.diagnostics.trex.diffraction_scan

or reduce the scan again without McStas::

    pixi run python -c "from mcdgs.diagnostics.trex.diffraction_scan import analyze; analyze()"
"""

from pathlib import Path

import numpy as np
import scipp as sc

from mcdgs.core.components import Single_crystal
from mcdgs.diagnostics import cli
from mcdgs.diagnostics import diffraction_scan as rotation_scan
from mcdgs.diagnostics.diffraction_scan import IN_PLANE, Q_STEP, QY_HALF_WIDTH
from mcdgs.diagnostics.trex import diffraction as plane
from mcdgs.diagnostics.trex import sample_runs as runs
from mcdgs.instruments.trex.model import DEFAULT_PARAMS, Params, TRex
from mcdgs.mcstas import Simulation, Workspace
from mcdgs.samples import la2ni7

NAME = "trex_diffraction_scan"
SCAN = f"{NAME}_omega_from_sample"  # the run folder of the scan, sample_pos -> detector
OMEGAS = np.arange(-30.0, 31.0, 2.0)  # deg, the angles of the scan
SCAN_SPLIT = 10  # McStas SPLIT on the crystal, at each angle
Q_LIMIT = (
    10.0  # 1/Å: the map spans -Q_LIMIT to Q_LIMIT in Qx and Qz, beyond 2 ki of every sub-pulse
)


# --8<-- [start:instrument]
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


def define_instrument() -> TRexWithRotatingLa2Ni7:
    """Return the instrument of this diagnostic: :class:`TRexWithRotatingLa2Ni7` at its defaults."""
    params = Params()  # the default control parameters
    return TRexWithRotatingLa2Ni7(params=params)


# --8<-- [end:instrument]


def simulate(omegas: np.ndarray = OMEGAS, beam_run: str = plane.TO_SAMPLE, seed: int = 1) -> Path:
    """Run the scan from ``sample_pos`` to the detector, once per angle, from a beam file.

    Parameters
    ----------
    omegas:
        The angles of the crystal, in deg.
    beam_run:
        The name of the run to the sample whose beam to start from: by default that of
        :mod:`.diagnostics.trex.diffraction`.
    seed:
        Seed of the random numbers; angle ``n`` takes ``seed + n * 1024``.

    Returns
    -------
    :
        The scan folder, ``SCAN`` in the runs folder of :meth:`.Workspace.find`, with one entry
        per angle.
    """
    sim = Simulation(define_instrument())
    scan = sim.scan_from_sample("omega", omegas, folder_name=SCAN, beam_run=beam_run, seed=seed)
    return scan.folder


def reduce_angles(
    scan: str | Path | None = None,
    to_sample: str | Path | None = None,
    step: float = Q_STEP,
    qy_half_width: float = QY_HALF_WIDTH,
) -> sc.DataArray:
    """Reduce every angle of a scan into its own (Qx, Qz) map, in the crystal's frame.

    The sub-pulses come from the run to the sample of :mod:`.diagnostics.trex.diffraction`, the same
    for every angle. Each angle is read and reduced on its own, so the scan never needs all its
    neutrons in memory at once.

    Parameters
    ----------
    scan:
        The scan folder of :func:`simulate`. If ``None``, ``SCAN`` in the runs folder of
        :meth:`.Workspace.find`.
    to_sample:
        The run to the sample. If ``None``, that of :mod:`.diagnostics.trex.diffraction`.
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
    workspace = Workspace.find()
    scan = Path(scan or workspace.run_folder(SCAN))
    to_sample = Path(to_sample or workspace.run_folder(plane.TO_SAMPLE))
    # the beam of the scan is that of diffraction, which its instrument ran
    pulses = runs.subpulse_timing_of_run(to_sample, plane.define_instrument())
    return rotation_scan.reduce_scan(
        scan, pulses, runs.read_detector, runs.pixel_solid_angle, Q_LIMIT, step, qy_half_width
    )


def analyze(
    out: Path | None = None,
    scan: str | Path | None = None,
    to_sample: str | Path | None = None,
) -> list[Path]:
    """Reduce a scan and save it as ``trex_diffraction_scan.gif`` and ``trex_diffraction_scan.png``.

    The GIF, from :func:`.diffraction_scan.animate`, adds one angle per frame; the figure, from
    :func:`.diffraction_scan.plot_scan`, is its last frame, the whole scan. It runs no McStas, so
    it can run again on the same scan.

    Parameters
    ----------
    out:
        Folder of the files. If ``None``, ``trex_diffraction_scan`` in the results folder of
        :meth:`.Workspace.find`.
    scan, to_sample:
        As :func:`reduce_angles`.

    Returns
    -------
    :
        The paths of the GIF and of the figure.
    """
    scan = Path(scan or Workspace.find().run_folder(SCAN))
    maps = rotation_scan.crop(reduce_angles(scan, to_sample))
    hkl = " ".join(str(i) for i in IN_PLANE)
    title = (
        f"T-Rex, La2Ni7 ({hkl}) in plane\n{define_instrument().running_conditions()}\n"
        "Q in the frame of the crystal at omega = 0"
    )
    folder = cli.output_folder(out, NAME)
    gif = rotation_scan.animate(maps, folder / f"{NAME}.gif", title)
    print(f"Saved {gif}")
    figure = rotation_scan.plot_scan(maps, title)
    return [gif, cli.save_figure(figure, folder / f"{NAME}.png")]


def main(out: Path | None = None) -> list[Path]:
    """Run the scan with :func:`simulate`, then reduce and draw it with :func:`analyze`.

    The beam file must exist: run :mod:`.diagnostics.trex.diffraction` first.

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
    main(out=cli.parse_out(__doc__))
