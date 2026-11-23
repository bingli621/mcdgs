import numpy as np
import pytest
import scipp as sc

from mcdgs.reduction.mdnorm import HBAR2_OVER_2M, momentum_transfer, q_e_coverage


def _coverage(
    two_theta: list[float], solid_angle: list[float], window: tuple[float, float]
) -> sc.Variable:
    """The coverage of pixels at 20 meV, on 0.5 meV bins of E and 0.05 1/Å bins of |Q|."""
    return q_e_coverage(
        sc.array(dims=["pixel"], values=two_theta, unit="deg"),
        sc.array(dims=["pixel"], values=solid_angle, unit="sr"),
        sc.scalar(20.0, unit="meV"),
        sc.linspace("energy_transfer", -20.0, 20.0, 81, unit="meV"),
        sc.linspace("Q", 0.0, 12.0, 241, unit="1/Å"),
        window,
    )


def test_hbar2_over_2m_is_that_of_the_neutron():
    assert HBAR2_OVER_2M == pytest.approx(2.0721, rel=1e-4)


def test_momentum_transfer_broadcasts_energy_against_angle():
    q = momentum_transfer(20.0, np.array([[0.0], [5.0]]), np.array([30.0, 60.0, 90.0]))
    assert q.shape == (2, 3)
    k = np.sqrt(20.0 / HBAR2_OVER_2M)
    assert q[0] == pytest.approx(2 * k * np.sin(np.radians([15.0, 30.0, 45.0])))


def test_each_bin_of_e_gets_the_solid_angle_times_its_length_in_the_window():
    coverage = _coverage([30.0, 120.0], [0.5, 0.25], window=(-0.8, 0.8))
    per_e = coverage.sum("Q").values  # sr meV in each bin of E
    inside = np.abs(sc.midpoints(sc.linspace("e", -20.0, 20.0, 81)).values) < 16.0
    assert per_e[inside] == pytest.approx(0.75 * 0.5)  # both pixels, the whole bin
    assert per_e[~inside] == pytest.approx(0.0)  # beyond ±0.8 Ei = ±16 meV
    assert coverage.sum().value == pytest.approx(0.75 * 32.0)


def test_the_elastic_bin_is_covered_at_the_q_of_its_angle():
    coverage = _coverage([90.0], [1.0], window=(-0.8, 0.8))
    elastic = coverage["energy_transfer", 40].values  # the bin from 0 to 0.5 meV
    k = np.sqrt(20.0 / HBAR2_OVER_2M)
    expected = momentum_transfer(20.0, 0.25, 90.0)  # at the middle of the bin
    assert np.flatnonzero(elastic)[0] <= int(expected / 0.05) <= np.flatnonzero(elastic)[-1]
    assert 2 * k * np.sin(np.radians(45.0)) == pytest.approx(momentum_transfer(20.0, 0.0, 90.0))


def test_uneven_bins_of_q_raise():
    with pytest.raises(ValueError, match="evenly spaced"):
        q_e_coverage(
            sc.array(dims=["pixel"], values=[90.0], unit="deg"),
            sc.array(dims=["pixel"], values=[1.0], unit="sr"),
            sc.scalar(20.0, unit="meV"),
            sc.array(dims=["energy_transfer"], values=[-1.0, 1.0], unit="meV"),
            sc.array(dims=["Q"], values=[0.0, 1.0, 3.0], unit="1/Å"),
            window=(-0.8, 0.8),
        )
