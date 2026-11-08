import pytest
import scipp as sc
from scippneutron.conversion.tof import energy_from_wavelength

from mcdgs.reduction.incident_energy import incident_flight

ORIGIN = sc.vector([0.0, 0.0, 0.0], unit="m")
SPEED_PER_ANGSTROM = 3956.034  # m/s times Å: the speed of a neutron of 1 Å


def _time(distance: float, wavelength: float) -> float:
    """The time in ms a neutron of ``wavelength`` Å takes over ``distance`` m."""
    return 1e3 * distance * wavelength / SPEED_PER_ANGSTROM


def test_a_flight_gives_its_length_time_wavelength_and_energy():
    end = sc.vector([0.0, 0.0, 3.3], unit="m")
    start_time, wavelength = 60.0, 1.8
    flight = incident_flight(
        ORIGIN,
        sc.scalar(start_time, unit="ms"),
        end,
        sc.scalar(start_time + _time(3.3, wavelength), unit="ms"),
    )
    assert flight["Ltotal"].to(unit="m").value == pytest.approx(3.3)
    assert flight["tof"].to(unit="ms").value == pytest.approx(_time(3.3, wavelength))
    assert flight["incident_wavelength"].to(unit="Å").value == pytest.approx(wavelength, rel=1e-6)
    expected = energy_from_wavelength(wavelength=sc.scalar(wavelength, unit="Å")).to(unit="meV")
    assert flight["incident_energy"].to(unit="meV").value == pytest.approx(expected.value, rel=1e-6)


def test_the_flight_path_is_the_straight_line_from_the_start_to_the_end():
    start = sc.vector([1.0, 0.0, 2.0], unit="m")
    end = sc.vector([4.0, 0.0, 6.0], unit="m")  # 3 m across, 4 m along: 5 m
    flight = incident_flight(start, sc.scalar(0.0, unit="ms"), end, sc.scalar(5.0, unit="ms"))
    assert flight["Ltotal"].to(unit="m").value == pytest.approx(5.0)


def test_one_start_time_holds_for_every_arrival():
    wavelengths = [1.5, 2.0, 2.5]
    arrivals = sc.array(
        dims=["subpulse"], values=[1.7 + _time(160.0, w) for w in wavelengths], unit="ms"
    )
    flight = incident_flight(
        ORIGIN, sc.scalar(1.7, unit="ms"), sc.vector([0.0, 0.0, 160.0], unit="m"), arrivals
    )
    assert flight["incident_energy"].dims == ("subpulse",)
    assert flight["incident_wavelength"].to(unit="Å").values == pytest.approx(wavelengths)


def test_a_start_time_per_neutron_and_one_end_time_broadcast_too():
    starts = sc.array(dims=["subpulse"], values=[0.0, 1.0], unit="ms")
    flight = incident_flight(
        ORIGIN, starts, sc.vector([0.0, 0.0, 10.0], unit="m"), sc.scalar(6.0, unit="ms")
    )
    assert flight["tof"].to(unit="ms").values == pytest.approx([6.0, 5.0])
