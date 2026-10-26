from dataclasses import dataclass

import pytest

from mcdgs.core.base import INSTRUMENTS, Instrument
from mcdgs.core.components import Arm, Progress_bar
from mcdgs.core.params import InstrParam, param


@dataclass(frozen=True)
class ToyParams:
    omega: float = param(0.0, "deg", "sample rotation")
    ei: float = param(25.0, "meV", "incident energy")


class Toy(Instrument):
    params_class = ToyParams
    name = "toy"

    def __init__(self, params: ToyParams | None = None) -> None:
        super().__init__(params)
        self.progress = Progress_bar()
        self.origin = Arm()
        self.sample = Arm(at=[0, 0, 1], rotated=[0, "omega", 0], relative="origin")
        self.detector = Arm(at=[0, 0, 2], relative="sample")

    @property
    def instr_params(self) -> list[InstrParam]:
        return [InstrParam("omega", self.params.omega, "deg", "sample rotation")]


def test_a_subclass_is_registered_by_its_class_name():
    assert INSTRUMENTS["Toy"] is Toy


def test_the_params_default_to_those_of_the_params_class():
    assert Toy().params == ToyParams()
    assert Toy(params=ToyParams(omega=5)).params.omega == 5


def test_update_params_replaces_the_frozen_params():
    toy = Toy()
    toy.update_params(omega=30)
    assert toy.params == ToyParams(omega=30)
    with pytest.raises(TypeError):
        toy.update_params(psi=1)


def test_components_are_named_after_their_attribute_in_beam_order():
    toy = Toy()
    assert list(toy.components) == ["progress", "origin", "sample", "detector"]
    assert toy.sample.name == "sample"


def test_an_attribute_set_to_none_is_no_component():
    toy = Toy()
    toy.progress = None
    assert "progress" not in toy.components


def test_a_component_cannot_take_a_reserved_name():
    toy = Toy()
    with pytest.raises(ValueError, match="positions: reserved name"):
        toy.positions = Arm()
    with pytest.raises(ValueError, match="_params: reserved name"):
        toy._params = Arm()


def test_positions_follow_the_instrument_parameters():
    toy = Toy()
    assert toy.positions()["detector"] == [0.0, 0.0, 3.0]
    toy.update_params(omega=90)
    assert toy.positions()["detector"] == [2.0, 0.0, 1.0]  # the sample's z is the lab's x


def test_instr_params_must_be_implemented():
    class Bare(Instrument):
        params_class = ToyParams
        name = "bare"

    with pytest.raises(NotImplementedError):
        _ = Bare().instr_params
    del INSTRUMENTS["Bare"]


def test_equality_compares_name_params_and_components():
    assert Toy() == Toy()
    assert Toy() != Toy(params=ToyParams(omega=1))
    other = Toy()
    other.detector = Arm(at=[0, 0, 4], relative="sample")
    assert Toy() != other
    renamed = Toy()
    renamed.name = "other"
    assert Toy() != renamed


def test_the_summary_shows_params_instr_params_and_components():
    text = repr(Toy())
    assert text.startswith("Toy 'toy'")
    assert "incident energy" in text
    assert "Instrument parameters (McStas):" in text
    assert "[0, omega, 0]" in text  # the sample's rotation, as given
