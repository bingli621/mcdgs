import pytest

from mcdgs.core.utils import format_value


@pytest.mark.parametrize(
    "value, text",
    [
        (1 / 3, "0.333333"),
        (2, "2"),
        (True, "True"),  # a bool is no number here
        ("omega", "omega"),
        (None, "None"),
        ([0, 0.5, "-chi"], "[0, 0.5, -chi]"),
    ],
)
def test_format_value(value, text):
    assert format_value(value) == text
