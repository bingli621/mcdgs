import pytest
import scipp as sc
from scippneutron.conversion.tof import energy_from_wavelength

from mcdgs.reduction.inelastic import QE_GRAPH, arrival_window, assign_subpulse, time_at_sample

ORIGIN = sc.vector([0.0, 0.0, 0.0], unit="m")
SPEED_PER_ANGSTROM = 3956.034  # m/s times Å: the speed of a neutron of 1 Å


def _time(distance: float, wavelength: float) -> float:
    """The time in ms a neutron of ``wavelength`` Å takes over ``distance`` m."""
    return 1e3 * distance * wavelength / SPEED_PER_ANGSTROM


def _energy(wavelength: float) -> sc.Variable:
    return energy_from_wavelength(wavelength=sc.scalar(wavelength, unit="Å")).to(unit="meV")


def _scatter(position: list[float], final: float, incident: float = 2.0) -> sc.DataArray:
    """One neutron of ``incident`` Å, scattered at the origin to ``position`` with ``final`` Å."""
    L2 = sc.norm(sc.vector(position, unit="m")).value
    return sc.DataArray(
        sc.scalar(1.0),
        coords={
            "source_position": sc.vector([0.0, 0.0, -160.0], unit="m"),
            "sample_position": ORIGIN,
            "position": sc.vector(position, unit="m"),
            "incident_energy": _energy(incident),
            "time_at_sample": sc.scalar(80.0, unit="ms"),
            "time_at_detector": sc.scalar(80.0 + _time(L2, final), unit="ms"),
        },
    )


def test_time_at_sample_adds_the_flight_from_the_monitor():
    t = time_at_sample(
        sc.scalar(50.0, unit="ms"),
        sc.vector([0.0, 0.0, -3.0], unit="m"),
        ORIGIN,
        _energy(1.5),
    )
    assert t.to(unit="ms").value == pytest.approx(50.0 + _time(3.0, 1.5), rel=1e-6)


def test_elastic_scattering_has_no_energy_transfer_and_q_from_the_angle():
    # 90° in the horizontal plane: |Q| = sqrt(2) k
    out = _scatter([3.0, 0.0, 0.0], final=2.0).transform_coords(
        ["energy_transfer", "Q", "Qx", "Qy", "Qz"], graph=QE_GRAPH
    )
    k = 2 * 3.141592653589793 / 2.0
    assert out.coords["energy_transfer"].to(unit="meV").value == pytest.approx(0.0, abs=1e-5)
    assert out.coords["Q"].to(unit="1/Å").value == pytest.approx(2**0.5 * k, rel=1e-6)
    # ki along z, kf along x: Q = ki - kf = (-k, 0, k)
    assert out.coords["Qx"].value == pytest.approx(-k, rel=1e-6)
    assert out.coords["Qy"].value == pytest.approx(0.0, abs=1e-9)
    assert out.coords["Qz"].value == pytest.approx(k, rel=1e-6)


def test_a_slower_neutron_has_lost_energy():
    out = _scatter([0.0, 0.0, 3.0], final=2.5).transform_coords(
        ["final_energy", "energy_transfer"], graph=QE_GRAPH
    )
    assert out.coords["final_energy"].to(unit="meV").value == pytest.approx(
        _energy(2.5).value, rel=1e-6
    )
    expected = (_energy(2.0) - _energy(2.5)).value
    assert out.coords["energy_transfer"].to(unit="meV").value == pytest.approx(expected, rel=1e-6)
    assert expected > 0


def test_the_arrival_window_spans_the_energy_transfer_range():
    ei = _energy(2.0)
    window = arrival_window(
        sc.scalar(10.0, unit="ms"), sc.scalar(3.0, unit="m"), ei, energy_transfer_ratio=(-3.0, 0.75)
    )
    # Ef = 4 Ei: half the wavelength; Ef = Ei / 4: twice the wavelength
    assert window["earliest_arrival"].value == pytest.approx(10.0 + _time(3.0, 1.0), rel=1e-6)
    assert window["latest_arrival"].value == pytest.approx(10.0 + _time(3.0, 4.0), rel=1e-6)


def test_the_arrival_window_broadcasts_sub_pulses_against_pixels():
    window = arrival_window(
        sc.array(dims=["subpulse"], values=[0.0, 5.0], unit="ms"),
        sc.array(dims=["pixel"], values=[3.0, 3.1, 3.2], unit="m"),
        sc.array(dims=["subpulse"], values=[_energy(2.0).value, _energy(1.8).value], unit="meV"),
    )
    assert set(window["earliest_arrival"].dims) == {"subpulse", "pixel"}


@pytest.mark.parametrize("ratio", [(0.5, 0.2), (-0.5, 1.0)])
def test_an_empty_or_unphysical_range_raises(ratio):
    with pytest.raises(ValueError, match="energy_transfer_ratio"):
        arrival_window(sc.scalar(0.0, unit="ms"), sc.scalar(3.0, unit="m"), _energy(2.0), ratio)


def _windows():
    # two pixels, with ids 7 and 3, and two sub-pulses: [10, 20) and [30, 40) ms, the second
    # pixel 1 ms later
    earliest = sc.array(dims=["pixel", "subpulse"], values=[[10.0, 30.0], [11.0, 31.0]], unit="ms")
    return sc.array(dims=["pixel"], values=[7, 3]), earliest, earliest + sc.scalar(10.0, unit="ms")


def test_each_neutron_gets_the_sub_pulse_of_its_window():
    ids, earliest, latest = _windows()
    out = assign_subpulse(
        sc.array(dims=["event"], values=[7, 7, 3, 3, 7, 3]),
        sc.array(dims=["event"], values=[15.0, 35.0, 10.5, 32.0, 25.0, 45.0], unit="ms"),
        ids,
        earliest,
        latest,
    )
    assert out["subpulse"].values.tolist() == [0, 1, 0, 1, 0, 1]
    # 10.5 ms is before the window of pixel 3, 25 ms between windows, 45 ms after the last
    assert out["inside"].values.tolist() == [True, True, False, True, False, False]


def test_the_times_may_be_in_another_unit():
    ids, earliest, latest = _windows()
    out = assign_subpulse(
        sc.array(dims=["event"], values=[3]),
        sc.array(dims=["event"], values=[0.035], unit="s"),
        ids,
        earliest,
        latest,
    )
    assert out["subpulse"].values.tolist() == [1]
    assert out["inside"].values.tolist() == [True]


def test_a_pixel_without_a_window_raises():
    ids, earliest, latest = _windows()
    with pytest.raises(ValueError, match="no arrival window"):
        assign_subpulse(
            sc.array(dims=["event"], values=[5]),
            sc.array(dims=["event"], values=[15.0], unit="ms"),
            ids,
            earliest,
            latest,
        )
