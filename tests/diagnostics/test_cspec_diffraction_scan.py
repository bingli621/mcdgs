from mcdgs.diagnostics.cspec.diffraction_scan import SCAN_SPLIT, define_instrument
from mcdgs.diagnostics.diffraction_scan import IN_PLANE


def test_the_crystal_turns_by_omega_on_the_instrument_of_the_scan():
    cspec = define_instrument()
    assert type(cspec.sample).__name__ == "Single_crystal"
    assert cspec.sample.split == SCAN_SPLIT
    assert list(cspec.sample.rotated) == [0, "omega", 0]
    assert "omega" in {param.name for param in cspec.instr_params}
    assert IN_PLANE == [2, 1, 0]


def test_a_histogram_monitor_follows_sample_pos_so_mcrun_counts_every_angle():
    names = list(define_instrument().components)
    after = names[names.index("sample_pos") + 1 : names.index("sample")]
    assert "sample_monitor" in after
