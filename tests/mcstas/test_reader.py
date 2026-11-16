import numpy as np
import pytest
import scipp as sc
from mcstasscript.data.data import McStasData, McStasMetaData

from mcdgs.mcstas import read_mcstas_map

BINS = {
    "toa": sc.linspace("toa", 0.0, 15.0, 4, unit="ms"),
    "wavelength": sc.linspace("wavelength", 0.0, 10.0, 3, unit="Å"),
}


def _monitor(name: str, shape: tuple[int, int]) -> McStasData:
    """A McStas monitor as McStasScript has it: (wavelength, time)."""
    metadata = McStasMetaData()
    metadata.component_name = name
    data = McStasData(metadata)
    data.Intensity = np.arange(np.prod(shape), dtype=float).reshape(shape)
    data.Error = np.ones(shape)
    return data


def test_the_map_of_the_named_monitor_is_read_as_time_by_wavelength():
    data = [_monitor("other", (2, 3)), _monitor("sample_monitor", (2, 3))]
    data[0].Intensity += 100
    intensity, error = read_mcstas_map("sample_monitor", data, BINS)
    assert intensity.shape == error.shape == (3, 2)  # time bins, wavelength bins
    np.testing.assert_array_equal(intensity, [[0, 3], [1, 4], [2, 5]])


def test_a_map_on_other_bins_is_refused():
    with pytest.raises(ValueError, match=r"sample_monitor: McStas has a map of \(2, 3\) bins"):
        read_mcstas_map("sample_monitor", [_monitor("sample_monitor", (3, 2))], BINS)


def test_a_name_of_several_data_sets_is_refused():
    data = [_monitor("sample_monitor", (2, 3)), _monitor("sample_monitor", (2, 3))]
    with pytest.raises(ValueError, match="sample_monitor: 2 data sets of that name"):
        read_mcstas_map("sample_monitor", data, BINS)
