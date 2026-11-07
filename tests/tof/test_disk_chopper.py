from mcdgs.core.components import NXdisk_chopper
from mcdgs.instruments.trex.model import TRex
from mcdgs.tof.builder import disk_chopper, get_choppers


def test_disk_takes_the_radius_and_the_slit_height():
    chopper = NXdisk_chopper(name="c", slit_edges=[0, 10, 90, 100], radius=0.3, yheight=0.05, nu=14)
    disk = disk_chopper(chopper, 10.0, {})
    assert disk.radius.value == 0.3
    assert disk.slit_height.values.tolist() == [0.05, 0.05]


def test_slits_without_a_height_reach_the_spindle():
    chopper = NXdisk_chopper(name="c", slit_edges=[0, 10], radius=0.3, nu=14)
    assert disk_chopper(chopper, 10.0, {}).slit_height.values.tolist() == [0.3]


def test_get_choppers_by_name_in_beam_order():
    choppers = get_choppers(TRex())
    assert list(choppers) == ["bw1", "bw2", "ps1", "ps2", "mc1", "mc2"]
