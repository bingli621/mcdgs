from dataclasses import dataclass

import pytest

from mcdgs.core.base import INSTRUMENTS, Instrument
from mcdgs.core.components import Arm, Guide_gravity, NXdisk_chopper, Progress_bar
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


def test_the_choppers_are_the_disk_choppers_in_beam_order():
    toy = Toy()
    toy.first = NXdisk_chopper()
    toy.second = NXdisk_chopper()
    assert toy.choppers == {"first": toy.first, "second": toy.second}
    assert toy.choppers["first"] is toy.first


def test_the_upstream_guide_of_a_chopper_is_the_last_guide_before_it():
    toy = Toy()
    toy.unguided = NXdisk_chopper()
    toy.guide_1 = Guide_gravity(w1=0.03, h1=0.1, l=1)
    toy.guide_2 = Guide_gravity(w1=0.03, h1=0.1, l=1)
    toy.arm = Arm()
    toy.first = NXdisk_chopper()
    toy.second = NXdisk_chopper()
    assert toy.upstream_guides == {"unguided": None, "first": toy.guide_2, "second": toy.guide_2}


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


def test_the_summary_aligns_values_to_the_left():
    lines = repr(Toy()).splitlines()
    header = lines[lines.index("Params:") + 1]
    omega, ei = lines[lines.index("Params:") + 3 : lines.index("Params:") + 5]
    assert omega[header.index("value")] == "0"
    assert ei[header.index("value")] == "2"


def test_the_summary_shows_params_instr_params_and_components():
    text = repr(Toy())
    assert text.startswith("Toy 'toy'")
    assert "incident energy" in text
    assert "Instrument parameters (McStas):" in text
    assert "[0, omega, 0]" in text  # the sample's rotation, as given


def test_chopper_frequencies_are_signed_by_name_in_beam_order():
    toy = Toy()
    toy.slow = NXdisk_chopper(slit_edges=[0, 10], nu=14, at=[0, 0, 3])
    toy.fast_1 = NXdisk_chopper(slit_edges=[0, 10], nu="omega", at=[0, 0, 4])
    toy.fast_2 = NXdisk_chopper(slit_edges=[0, 10], nu="-2*omega", at=[0, 0, 5])
    toy.update_params(omega=70)
    assert toy.chopper_frequencies() == {"slow": 14.0, "fast_1": 70.0, "fast_2": -140.0}
    assert toy.chopper_frequencies("fast") == {"fast_1": 70.0, "fast_2": -140.0}
    assert toy.chopper_frequencies("none") == {}


def test_running_conditions_default_to_the_name():
    assert Toy().running_conditions() == "toy"


def test_numeric_instr_params_leave_out_strings():
    class Named(Toy):
        @property
        def instr_params(self) -> list[InstrParam]:
            return [*super().instr_params, InstrParam("label", "a", "", "a string")]

    assert Named(params=ToyParams(omega=5)).numeric_instr_params() == {"omega": 5}
    del INSTRUMENTS["Named"]
