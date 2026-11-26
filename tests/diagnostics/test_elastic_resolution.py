import numpy as np
import pytest
import scipp as sc

from mcdgs.diagnostics.elastic_resolution import scattered_fractions


def test_each_neutron_is_its_fraction_of_the_beam_divided_by_kf_over_ki():
    events = sc.DataArray(
        sc.array(dims=["event"], values=[2.0, 2.0, 2.0], unit="1/sr"),
        coords={
            "subpulse": sc.array(dims=["event"], values=[0, 1, 1]),
            "energy_transfer": sc.array(dims=["event"], values=[0.0, 7.5, -30.0], unit="meV"),
            "solid_angle": sc.array(dims=["event"], values=[0.5, 0.5, 0.5], unit="sr"),
        },
    )
    pulses = sc.DataArray(
        sc.ones(dims=["subpulse"], shape=[2]),
        coords={"incident_energy": sc.array(dims=["subpulse"], values=[5.0, 10.0], unit="meV")},
    )
    fractions = scattered_fractions(events, pulses)
    assert fractions.unit == sc.units.one
    # elastic: as it is; kf/ki = 1/2 on the loss side, 2 on the gain side, with the Ei of 10 meV
    assert fractions.values == pytest.approx(np.array([1.0, 2.0, 0.5]))
