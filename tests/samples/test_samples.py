from dataclasses import fields
from pathlib import Path

import numpy as np
import pytest

from mcdgs.core.components import Incoherent, Phonon_simple, Single_crystal
from mcdgs.core.components.crystal import CIF_DIR
from mcdgs.samples import aluminum, la2ni7, lead, vanadium

ALUMINUM = {  # the parameters of aluminum()
    "xwidth": 0.005, "yheight": 0.005, "zdepth": 0.005, "order": 1, "mosaic": 5.0,
    "reflections": "Al.lau", "lattice_params": [4.0495] * 3 + [90.0] * 3,
    "along_beam": [0, 0, 1], "in_plane": [1, 0, 0],
}  # fmt: skip


def test_la2ni7_orientation():
    s = la2ni7()
    a, b, c = (np.array(v) for v in (s.astar, s.bstar, s.cstar))
    # (001) along +z, the beam
    assert c[:2] == pytest.approx(0, abs=1e-9) and c[2] > 0
    # (1.357, 1, 0) in the horizontal plane (y = 0), toward +x
    q = 1.357 * a + b
    assert q[1] == pytest.approx(0, abs=1e-9) and q[0] > 0
    # right-handed lattice
    assert np.dot(np.cross(a, b), c) > 0


def test_the_samples_are_crystals_placed_like_any_component():
    for make in (aluminum, la2ni7):
        crystal = make(name="sample", at=[0, 0, 0.1], rotated=[0, 5, 0], relative="arm")
        assert type(crystal) is Single_crystal
        assert crystal.name == "sample"
        assert (crystal.at, crystal.rotated, crystal.relative) == ([0, 0, 0.1], [0, 5, 0], "arm")
        assert (make().at, make().relative) == ([0, 0, 0], None)
    assert la2ni7(mosaic=60).mosaic == 60


def test_the_cif_of_a_crystal_is_a_path_inside_the_package():
    assert not Path(la2ni7().reflections).is_absolute()
    assert (CIF_DIR / la2ni7().reflections).is_file()


def test_the_mcstas_parameters_give_the_cif_as_a_full_path():
    assert la2ni7().mcstas_parameters()["reflections"] == str(CIF_DIR / "La2Ni7.cif")
    assert aluminum().mcstas_parameters()["reflections"] == "Al.lau"  # McStas's own


def test_the_mcstas_parameters_give_the_reciprocal_vectors_in_the_lab_frame():
    crystal = aluminum()
    parameters = crystal.mcstas_parameters()
    assert parameters["recip_cell"] == 1
    for axis, vector in zip("abc", (crystal.astar, crystal.bstar, crystal.cstar)):
        assert [parameters[f"{axis}{x}"] for x in "xyz"] == vector
    assert parameters["cz"] == pytest.approx(2 * np.pi / 4.0495)


def test_the_crystal_keeps_its_lattice_params():
    crystal = la2ni7()
    assert crystal.lattice_params == [5.0556, 5.0556, 24.598, 90.0, 90.0, 120.0]
    B = crystal.reciprocal_basis()  # columns a*, b*, c* in the crystal frame
    assert np.linalg.norm(B[:, 2]) == pytest.approx(np.linalg.norm(crystal.cstar))


def test_the_reciprocal_vectors_are_not_parameters():
    crystal = Single_crystal(**ALUMINUM)
    assert crystal == aluminum()
    assert "astar" not in [f.name for f in fields(crystal)]  # derived, so never in a file


def test_the_reciprocal_vectors_follow_a_change_of_the_crystal():
    crystal = aluminum()
    before = crystal.astar
    crystal.along_beam, crystal.in_plane = [1, 0, 0], [0, 1, 0]
    assert crystal.astar != before
    assert crystal.astar == pytest.approx([0, 0, 2 * np.pi / 4.0495])  # a* now along the beam


def test_size_mosaic_and_order_default_to_mcstas():
    crystal = Single_crystal(
        reflections="Al.lau", lattice_params=[4.0495] * 3 + [90.0] * 3,
        along_beam=[0, 0, 1], in_plane=[1, 0, 0],
    )  # fmt: skip
    assert (crystal.xwidth, crystal.yheight, crystal.zdepth, crystal.mosaic, crystal.order) == (
        0, 0, 0, -1, 0,
    )  # fmt: skip


def test_lead_is_a_phonon_sample_placed_like_any_component():
    sample = lead(name="sample", at=[0, 0, 0.1], rotated=[0, 5, 0], relative="arm")
    assert type(sample) is Phonon_simple
    assert sample.name == "sample"
    assert (sample.at, sample.rotated, sample.relative) == ([0, 0, 0.1], [0, 5, 0], "arm")
    assert (lead().at, lead().relative) == ([0, 0, 0], None)


def test_lead_has_the_lattice_the_atom_and_the_longitudinal_sound_velocity_of_lead():
    sample = lead()
    assert (sample.a, sample.b, sample.M) == (4.9508, 9.405, 207.2)
    assert (sample.sigma_abs, sample.sigma_inc) == (0.171, 0.0030)
    assert sample.c == pytest.approx(13.75, abs=0.005)  # 2089 m/s, in meV Å
    assert sample.branch_maximum == pytest.approx(11.11, abs=0.005)  # meV


def test_the_focusing_window_of_lead_is_given_by_the_caller():
    assert (lead().focus_r, lead().focus_aw, lead().focus_ah, lead().target_index) == (0, 0, 0, 0)
    sample = lead(target_index=2, focus_aw=60, focus_ah=10, c=7.54, T=100)
    assert (sample.target_index, sample.focus_aw, sample.focus_ah) == (2, 60, 10)
    assert (sample.c, sample.T) == (7.54, 100)


def test_vanadium_is_an_incoherent_sample_placed_like_any_component():
    sample = vanadium(name="sample", at=[0, 0, 0.1], rotated=[0, 5, 0], relative="arm")
    assert type(sample) is Incoherent
    assert sample.name == "sample"
    assert (sample.at, sample.rotated, sample.relative) == ([0, 0, 0.1], [0, 5, 0], "arm")
    assert (vanadium().at, vanadium().relative) == ([0, 0, 0], None)


def test_vanadium_is_a_solid_cylinder_with_the_cross_sections_of_vanadium():
    sample = vanadium()
    assert (sample.radius, sample.yheight, sample.thickness, sample.pack) == (0.005, 0.01, 0, 1)
    assert (sample.sigma_abs, sample.sigma_inc, sample.Vc) == (5.08, 5.08, 13.827)
    assert sample.Vc == pytest.approx(3.024**3 / 2, abs=5e-4)  # two atoms per bcc cell
    assert sample.scattering_length == pytest.approx(0.0272, abs=5e-5)  # m
    assert (sample.f_QE, sample.deltaE) == (0, 0)  # elastic


def test_the_focusing_window_of_vanadium_is_given_by_the_caller():
    sample = vanadium()
    assert (sample.focus_r, sample.focus_aw, sample.focus_ah, sample.target_index) == (0, 0, 0, 0)
    sample = vanadium(target_index=2, focus_aw=60, focus_ah=10, order=1, pack=0.6)
    assert (sample.target_index, sample.focus_aw, sample.focus_ah) == (2, 60, 10)
    assert (sample.order, sample.pack) == (1, 0.6)
