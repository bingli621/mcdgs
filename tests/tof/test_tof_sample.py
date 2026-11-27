from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
import scipp as sc
import tof
from scipp.constants import h, m_n

from mcdgs.core.components import Guide_gravity, Incoherent, Isotropic_Sqw
from mcdgs.instruments.trex.model import TRex
from mcdgs.samples import al2o3, lead, magnon_powder, vanadium, write_magnon_sqw
from mcdgs.samples.inelastic.powder import MAGNON_GAP, MAGNON_TOP
from mcdgs.tof import beamline, build_tof
from mcdgs.tof.builder import DETECTOR
from mcdgs.tof.sample import (
    QUASIELASTIC_GAMMAS,
    energy_transfer_band,
    inelastic_sample,
    sqw_energy_transfers,
    two_band_final_energy,
)


def _write_sqw(path: Path, w: list[float], sqw: list[list[float]]) -> Path:
    q = [0.0, 1.0][: len(sqw)]
    rows = "\n".join(" ".join(str(v) for v in row) for row in sqw)
    path.write_text(f"# q\n{' '.join(map(str, q))}\n# w\n{' '.join(map(str, w))}\n# sqw\n{rows}\n")
    return path


@pytest.fixture(scope="module")
def magnon_file(tmp_path_factory) -> Path:
    # |Q| up to 2 1/Å, past the corner of the zone, finely enough to reach its top
    return write_magnon_sqw(tmp_path_factory.mktemp("magnon") / "m.sqw", q_max=2.0)


def test_the_band_is_where_sqw_reaches_the_threshold(tmp_path):
    path = _write_sqw(
        tmp_path / "s.sqw", [-3, -2, -1, 0, 1, 2, 3], [[0, 1e-30, 1, 0, 1, 1, 0], [0] * 7]
    )
    assert sqw_energy_transfers(path, threshold=1e-20) == (1.0, 2.0)
    assert sqw_energy_transfers(path, threshold=0) == (1.0, 2.0)  # 1e-30 at w=-2: within it
    assert sqw_energy_transfers(path, threshold=2) is None


def test_a_file_of_one_side_gives_the_same_band(tmp_path):
    path = _write_sqw(tmp_path / "s.sqw", [0, 1, 2, 3], [[0, 0, 1, 1]])
    assert sqw_energy_transfers(path, threshold=1e-20) == (2.0, 3.0)


def test_a_powder_file_is_elastic(tmp_path):
    path = tmp_path / "p.laz"
    path.write_text("# lines\n1 2 3\n4 5 6\n")
    assert sqw_energy_transfers(path, threshold=1e-20) == (0.0, 0.0)


def test_a_file_of_two_blocks_raises(tmp_path):
    path = tmp_path / "s.sqw"
    path.write_text("# q\n0 1\n# w\n0 1\n")
    with pytest.raises(ValueError, match="2 blocks"):
        sqw_energy_transfers(path, threshold=1e-20)


def test_the_magnon_band_is_its_gap_to_its_top_widened_by_its_line(magnon_file):
    w_min, w_max = energy_transfer_band(magnon_powder(magnon_file))
    assert MAGNON_GAP - 0.2 < w_min < MAGNON_GAP
    assert MAGNON_TOP < w_max < MAGNON_TOP + 0.2


def test_the_sample_takes_the_widest_band_of_its_parts(tmp_path):
    coh = _write_sqw(tmp_path / "c.sqw", [0, 1, 2, 3], [[0, 0, 1, 1]])
    inc = _write_sqw(tmp_path / "i.sqw", [0, 1, 2, 3], [[0, 1, 1, 0]])
    sample = Isotropic_Sqw(Sqw_coh=str(coh), Sqw_inc=str(inc), sigma_coh=1, sigma_inc=1)
    assert energy_transfer_band(sample) == (1.0, 3.0)
    assert energy_transfer_band(replace(sample, sigma_inc=-1)) == (2.0, 3.0)  # left out
    # the incoherent part with no file is elastic
    assert energy_transfer_band(replace(sample, Sqw_inc=None)) == (0.0, 3.0)
    # elastic only: no band
    assert energy_transfer_band(replace(sample, Sqw_coh=None, Sqw_inc=None)) is None


def test_the_phonon_band_is_its_branch_from_the_zone_center_to_the_boundary():
    sample = lead()
    assert energy_transfer_band(sample) == (0.0, pytest.approx(sample.branch_maximum))
    w_min, w_max = energy_transfer_band(replace(sample, gap=3.0))
    assert (w_min, w_max) == (3.0, pytest.approx(np.hypot(3.0, sample.branch_maximum)))


@pytest.mark.parametrize(
    ("params", "band"),
    [
        ({"Etrans": 10, "deltaE": 2}, (8.0, 12.0)),
        ({"Etrans": -10, "deltaE": -2}, (8.0, 12.0)),  # gain, and a negative width
        ({"Etrans": 1, "deltaE": 2}, (0.0, 3.0)),  # across the elastic line
        ({"f_QE": 0.5, "gamma": 0.2}, (0.0, QUASIELASTIC_GAMMAS * 0.2)),
        ({"Etrans": 10, "deltaE": 2, "f_QE": 0.5, "gamma": 0.2}, (8.0, 12.0)),  # QE unused
    ],
)
def test_the_incoherent_band_follows_mcstas(params, band):
    assert energy_transfer_band(Incoherent(**params)) == pytest.approx(band)


@pytest.mark.parametrize(
    "component",
    [vanadium(), Incoherent(f_QE=0.5), al2o3(), Guide_gravity(w1=0.05, h1=0.05, l=1)],
    ids=["vanadium", "QE without width", "powder", "guide"],
)
def test_elastic_samples_and_other_components_get_no_tof_sample(component):
    assert energy_transfer_band(component) is None
    assert inelastic_sample("sample", component, 160.0, seed=1) is None


def test_the_transfers_fill_both_bands_as_often():
    final = two_band_final_energy(5.0, 20.0, np.random.default_rng(1))
    ei = sc.full(dims=["event"], shape=[100_000], value=50.0, unit="meV")
    w = (ei - final(ei)).values
    assert ((np.abs(w) >= 5) & (np.abs(w) <= 20)).all()
    assert (w > 0).mean() == pytest.approx(0.5, abs=0.01)
    assert np.abs(w).mean() == pytest.approx(12.5, abs=0.1)  # uniform in |w|


def test_trex_with_the_magnon_gets_a_tof_sample_at_the_sample(magnon_file):
    trex = TRex(magnon_powder(magnon_file))
    model = build_tof(trex, neutrons=1000, seed=1)
    sample = model.components["sample"]
    assert isinstance(sample, tof.InelasticSample)
    distance = {n: d for n, _, d in beamline(trex)[1]}["sample_pos"]
    assert sample.distance.value == pytest.approx(distance, abs=1e-12)


@pytest.mark.parametrize(
    "sample", [lead(), Incoherent(Etrans=5, deltaE=1), vanadium()], ids=["lead", "Etrans", "V"]
)
def test_trex_gets_a_tof_sample_for_an_inelastic_sample_only(sample):
    model = build_tof(TRex(sample), neutrons=1000, seed=1)
    inelastic = energy_transfer_band(sample) is not None
    assert isinstance(model.components.get("sample"), tof.InelasticSample) is inelastic


def _energy(wavelength: sc.Variable) -> np.ndarray:
    return (h**2 / (2 * m_n * wavelength**2)).to(unit="meV").values


def test_each_neutron_reaches_the_detector_with_a_transfer_within_the_bands(magnon_file):
    trex = TRex(magnon_powder(magnon_file))
    w_min, w_max = energy_transfer_band(trex.sample)
    result = build_tof(trex, neutrons=20_000, seed=1).run()
    before = result["sample_monitor"].data["pulse", 0]
    after = result[DETECTOR].data["pulse", 0]
    # through the choppers, and left with an energy above 0
    reached = ~after.masks["blocked_by_others"].values & ~np.isnan(
        after.coords["wavelength"].values
    )
    w = _energy(before.coords["wavelength"])[reached] - _energy(after.coords["wavelength"])[reached]
    assert reached.sum() > 20  # the monochromatic choppers let few through
    assert ((np.abs(w) >= w_min - 1e-9) & (np.abs(w) <= w_max + 1e-9)).all()
    assert (w > 0).any() and (w < 0).any()


def test_the_same_seed_gives_the_same_transfers(magnon_file):
    trex = TRex(magnon_powder(magnon_file))
    first, second = (build_tof(trex, neutrons=1000, seed=3).run() for _ in range(2))
    # tof gives a neutron left with an energy below 0 a wavelength of NaN
    assert sc.identical(first[DETECTOR].data, second[DETECTOR].data, equal_nan=True)
