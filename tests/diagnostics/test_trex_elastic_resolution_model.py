from dataclasses import replace

import numpy as np
import pytest
import scipp as sc

from mcdgs.diagnostics.trex.elastic_resolution_incoherent_short import (
    MIN_TWO_THETA,
    TRexElasticResolution,
    sample,
)
from mcdgs.diagnostics.trex.elastic_resolution_model import (
    ChopperPair,
    elastic_resolution,
    energy_of,
    fwhm_of,
    geometry_of,
    literature_resolution,
    opening_of,
    optimal_ratio,
    pair_transmission,
    read_pair,
    sample_path_density,
    slit_opening,
)
from mcdgs.instruments.trex.detectors import BANANA_RADIUS, BANANA_THETA_LIMITS
from mcdgs.instruments.trex.model import Params, TRex

TWO_THETA = ((BANANA_THETA_LIMITS[0], -MIN_TWO_THETA), (MIN_TWO_THETA, BANANA_THETA_LIMITS[1]))


def geometry(mode):
    # at the default lambda0, 2 Å, of the McStas widths below, not the 1 Å of the diagnostic
    params = Params(mode=mode)
    radius = float(sample().radius)
    return geometry_of(TRexElasticResolution(params=params), BANANA_RADIUS, radius, TWO_THETA)


def half_max_width(times, values):
    above = times[values >= values.max() / 2]
    return above[-1] - above[0]


def test_a_pair_over_a_beam_wider_than_its_slit_is_open_as_a_triangle_of_half_the_slit_time():
    # eq. (1) of Gaspar (2007): the FWHM of a counter-rotating pair is a / (2 u)
    pair = ChopperPair(distance=0, spacing=0, slit_arc=0.02, rim_speed=400, beam_width=0.03)
    times = np.linspace(-1e-4, 1e-4, 200001)
    opening = pair_transmission(pair, 0.0, times)
    assert half_max_width(times, opening) == pytest.approx(0.02 / (2 * 400), rel=1e-3)
    assert opening[np.abs(times) >= 0.02 / (2 * 400)].max() == 0  # the slits no longer overlap


def test_a_beam_narrower_than_the_slit_opens_the_pair_longer_up_to_one_disk():
    wide = ChopperPair(distance=0, spacing=0, slit_arc=0.0153, rim_speed=554, beam_width=0.0224)
    times = np.linspace(-1e-4, 1e-4, 200001)
    widths = [
        half_max_width(times, pair_transmission(replace(wide, beam_width=beam), 0.0, times))
        for beam in (0.0224, 0.0153, 0.005, 1e-6)
    ]
    assert widths == sorted(widths)  # longer as the beam narrows
    assert widths[0] == pytest.approx(0.0153 / (2 * 554), rel=1e-3)  # the triangle
    assert widths[-1] == pytest.approx(0.0153 / 554, rel=1e-3)  # a single disk


def test_the_offset_between_the_disks_moves_the_opening_but_not_its_area():
    # on a beam much wider than the slits, the window where they overlap stays on it
    pair = ChopperPair(distance=0, spacing=0.1, slit_arc=0.12, rim_speed=416, beam_width=1.0)
    times = np.linspace(-1e-3, 1e-3, 200001)
    centered = pair_transmission(pair, 0.0, times)
    shifted = pair_transmission(pair, 2e-5, times)
    assert np.average(times, weights=shifted) == pytest.approx(1e-5, abs=1e-7)
    assert shifted.sum() == pytest.approx(centered.sum(), rel=1e-3)


def test_the_sample_term_is_a_semicircle_for_one_angle():
    paths = np.linspace(-0.02, 0.02, 4001)
    density = sample_path_density(0.01, ((89.99, 90.01),), paths)
    reach = 2 * 0.01 * np.sin(np.deg2rad(45))  # 2 r sin θ
    assert density[np.abs(paths) > reach * 1.001].max() == 0
    assert density[2000] == pytest.approx(2 / (np.pi * reach), rel=1e-3)  # its height at 0


def test_the_fwhm_of_a_gaussian_is_its_fwhm():
    x = np.linspace(-3, 3, 2001)
    assert fwhm_of(x, np.exp(-(x**2) / (2 * 0.4**2))) == pytest.approx(2.3548 * 0.4, rel=1e-4)


def test_the_flight_time_weighs_the_monochromatic_offset_one_more_than_the_pulse_shaping():
    a, b = geometry("FF").coefficients
    assert a - b == pytest.approx(1)
    assert b == pytest.approx((1.767 + BANANA_RADIUS) / 54.03, rel=1e-3)


def test_the_pairs_of_the_mode_take_its_slits():
    ff, rr = geometry("FF"), geometry("RR")
    assert ff.pulse_shaping.slit_arc == pytest.approx(np.deg2rad(35) * 0.35)
    assert rr.pulse_shaping.slit_arc == pytest.approx(np.deg2rad(20) * 0.35)
    assert ff.monochromatic.slit_arc == pytest.approx(np.deg2rad(4.4) * 0.35)
    assert rr.monochromatic.slit_arc == pytest.approx(np.deg2rad(2.5) * 0.35)
    assert ff.monochromatic.beam_width == pytest.approx(0.022433, abs=1e-6)


def test_parked_choppers_have_no_resolution():
    params = Params(rrm=0)
    with pytest.raises(ValueError, match="parked"):
        read_pair(TRex(params=params, source_masked=False), ("ps1", "ps2"))


@pytest.mark.parametrize(
    ("mode", "fwhm"),
    [  # the McStas FWHM at 20.45 meV and 15.59 meV, of the vanadium at the default 2 Å
        ("FF", (1.061, 0.720)),
        ("FR", (0.954, 0.632)),
        ("RF", (0.846, 0.571)),
        ("RR", (0.592, 0.402)),
    ],
)
def test_the_model_gives_the_width_of_the_simulated_lines_within_five_percent(mode, fwhm):
    energies = sc.array(dims=["subpulse"], values=[20.45, 15.59], unit="meV")
    widths = elastic_resolution(geometry(mode), energies)["fwhm"].values
    assert widths == pytest.approx(fwhm, rel=0.05)


def test_the_narrow_pulse_shaping_slits_do_more_than_the_narrow_monochromatic_ones():
    energies = sc.array(dims=["subpulse"], values=[20.45], unit="meV")
    widths = {mode: elastic_resolution(geometry(mode), energies) for mode in ("FF", "FR", "RF")}
    assert widths["RF"]["fwhm"].values[0] < widths["FR"]["fwhm"].values[0]
    assert widths["FF"]["pulse_shaping"].values[0] > widths["FF"]["monochromatic"].values[0]


def test_the_formula_of_lechner_adds_its_three_terms_in_quadrature():
    energies = sc.array(dims=["subpulse"], values=[20.45], unit="meV")
    widths = literature_resolution(geometry("FF"), energies)
    terms = [widths[name].values[0] for name in ("pulse_shaping", "monochromatic", "sample")]
    assert widths["fwhm"].values[0] == pytest.approx(np.sqrt(np.sum(np.square(terms))))


def test_the_formula_of_lechner_has_the_prefactor_of_his_eq_3():
    # 647.2 / λ³ µeV·m/µs is 2 Ei v0: at 2 Å, τ2 = 10 µs on 1 m to the detector and L_PM = 1 m
    g = geometry("FF")
    pair = replace(g.monochromatic, slit_arc=2 * g.monochromatic.rim_speed * 10e-6)  # τ2 = 10 µs
    tiny = replace(g.pulse_shaping, slit_arc=1e-12)
    simple = replace(
        g,
        pulse_shaping=replace(tiny, distance=0.0),
        monochromatic=replace(pair, distance=1.0),
        sample_distance=0.0,
        detector_distance=1.0,
        sample_radius=1e-12,
    )
    energies = sc.array(dims=["subpulse"], values=[energy_of(2.0)], unit="meV")
    width = literature_resolution(simple, energies)["fwhm"].values[0]
    # eq. (3): 647.2 · τ2 (L12 + LSD) / (L12 LSD λ³) µeV
    assert width * 1000 == pytest.approx(647.2 * 10 * 2 / 2**3, rel=1e-3)


def test_the_slit_opening_is_the_triangle_of_the_literature():
    g = geometry("RR")
    assert slit_opening(g.monochromatic) == pytest.approx(opening_of(g.monochromatic)[0], rel=1e-3)
    assert slit_opening(g.pulse_shaping) < opening_of(g.pulse_shaping)[0]  # a narrower beam


def test_the_optimal_ratio_of_lechner_makes_the_two_chopper_terms_equal():
    g = geometry("FF")
    a, b = g.coefficients
    assert optimal_ratio(g) == pytest.approx(a / b)
    assert optimal_ratio(g) == pytest.approx(12.33, abs=0.01)


@pytest.mark.parametrize(("mode", "flux"), [("FR", 0.33), ("RF", 0.50), ("RR", 0.16)])
def test_the_areas_of_the_openings_give_the_flux_of_mcstas_over_ff(mode, flux):
    # the McStas flux at the pre-sample monitor over that of FF, at 20.45 meV
    def areas(m):
        g = geometry(m)
        return opening_of(g.pulse_shaping)[1] * opening_of(g.monochromatic)[1]

    assert areas(mode) / areas("FF") == pytest.approx(flux, abs=0.015)
