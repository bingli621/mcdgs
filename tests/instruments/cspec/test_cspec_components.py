import numpy as np
import pytest

from mcdgs.core import place_absolute
from mcdgs.core.components import Component, Guide_gravity
from mcdgs.instruments.cspec.choppers import DISKS, chopper
from mcdgs.instruments.cspec.guides.bender import guides_bender
from mcdgs.instruments.cspec.guides.frame import arm_iscs
from mcdgs.instruments.cspec.guides.straight import SAMPLE_PATH, guides_straight, sample_pos
from mcdgs.instruments.cspec.sources import source, wavelength_band


def _components() -> dict[str, Component]:
    """The guides, gaps, choppers and sample position, in beam order, without the source."""
    components: dict[str, Component] = {"arm_iscs": arm_iscs()}
    for piece in guides_bender() + guides_straight():
        components |= piece
    components["sample_pos"] = sample_pos()
    for name in DISKS:
        components[name], _ = chopper(name, None)
    return components


@pytest.fixture(scope="module")
def placed():
    return place_absolute(_components(), {})


def _exit(placed, components, name):
    position, rotation = placed[name]
    return position + rotation.T @ np.array([0, 0, components[name].l])


# the McStas start and end of a few elements, cells N to S of the Specification sheet
SHEET_MCSTAS = {
    "guide_01_02": (  # row 17
        (0.007582830807791673, -7.828132780329196e-05, 2.9042156263943757),
        (-0.0008852356922715923, -0.0007026562871538789, 3.9031795406603944),
    ),
    "guide_06_012": (  # row 52
        (-0.13036649621519725, -0.21117549160945825, 28.402167260328625),
        (-0.1342865532859437, -0.24215412889178342, 30.38392128161384),
    ),
    "guide_06_14": (  # row 78
        (-0.10374510652956681, -0.4224208564451449, 54.40196565200603),
        (-0.09519842506217713, -0.4224208564451449, 56.40194739048163),
    ),
    "guide_08_01": (  # row 135
        (0.1157262281918977, -0.4224208564451449, 105.75975782898323),
        (0.126324967879608, -0.4224208564451449, 108.23993518286682),
    ),
    "guide_09_02": (  # row 189
        (0.34121106071431906, -0.4224208564451449, 158.52474601265263),
        (0.3424930629344276, -0.4224208564451449, 158.824743273424),
    ),
}


@pytest.mark.parametrize("name", SHEET_MCSTAS)
def test_the_guides_start_and_end_where_the_sheet_has_them_in_mcstas(placed, name):
    # the sheet's McStas columns round the ISCS offset to 0.1 mm: 22 µm
    start, end = SHEET_MCSTAS[name]
    components = _components()
    np.testing.assert_allclose(placed[name][0], start, atol=3e-5)
    np.testing.assert_allclose(_exit(placed, components, name), end, atol=3e-5)


def test_the_guide_system_starts_at_its_onset_in_mcstas(placed):
    # cells B7, C7, D7 of Definitions
    np.testing.assert_allclose(placed["guide_01_011"][0], (0.0163, 0, 1.90424), atol=1e-12)


def test_each_guide_ends_where_the_next_element_starts():
    # the gaps between the elements are 1 to 40 mm, and the arm of a gap starts at the guide exit
    components = _components()
    placed = place_absolute(components, {})
    guides = [n for n, c in components.items() if isinstance(c, Guide_gravity)]
    starts = [n for n in components if n.startswith(("guide_", "arm_"))]
    for name in guides:
        following = starts[starts.index(name) + 1 :]
        if not following:
            continue
        gap = np.linalg.norm(_exit(placed, components, name) - placed[following[0]][0])
        assert gap < 0.0401, name  # the longest, row 50, is 40 mm


# the center of the first disk of each chopper, in ISCS, cells B to D of Definitions, in mm
DEFINITIONS_ISCS = {
    "bw1": (14908.65, 19.62, -49.059),
    "bw2": (20430.39, 40.73, -101.825),
    "bw3": (104579.53, 990.76, -422.42),
    "ps1": (105636.03, 1004.585, -422.42),
    "mc1": (158466.88437698322, 1691.145, -422.42),
}


@pytest.mark.parametrize("name", DEFINITIONS_ISCS)
def test_the_first_disk_of_each_chopper_is_where_the_definitions_have_it(placed, name):
    x, y, z = DEFINITIONS_ISCS[name]
    origin, rotation = placed["arm_iscs"]
    expected = origin + rotation.T @ (np.array([y, z, x]) / 1000)
    # the definitions round the lateral position: m1y is 0.2 mm off the trajectory
    np.testing.assert_allclose(placed[name][0], expected, atol=2.5e-4)


def test_the_sample_is_on_the_trajectory_of_the_last_guide(placed):
    position, rotation = placed["guide_09_02"]
    offset = rotation @ (placed["sample_pos"][0] - position)
    np.testing.assert_allclose(offset[:2], 0, atol=2.5e-4)
    np.testing.assert_allclose(placed["sample_pos"][1], rotation, atol=1e-12)


def test_the_sample_is_160_m_from_the_moderator():
    assert SAMPLE_PATH == pytest.approx(160.004, abs=1e-3)


def test_the_source_focuses_on_the_whole_entrance_of_the_guide_system():
    s = source()
    entrance = _components()["guide_01_011"]
    assert s.dist == pytest.approx(1.90424)
    assert s.focus_xw / 2 == pytest.approx(0.0163 + entrance.w1 / 2)
    assert s.focus_yh == entrance.h1


def test_the_band_of_the_source_is_centered_on_lambda0():
    l_min, l_max = wavelength_band(4.0)
    assert (l_min + l_max) / 2 == pytest.approx(4.0)
    assert l_max - l_min == pytest.approx(1.853, abs=1e-3)


def test_the_second_disk_of_each_pair_is_10_mm_behind_the_first():
    assert DISKS["ps2"].distance - DISKS["ps1"].distance == pytest.approx(0.010)
    assert DISKS["mc2"].distance - DISKS["mc1"].distance == pytest.approx(0.010)
