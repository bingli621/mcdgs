from dataclasses import dataclass

import pytest

from mcdgs.core.components import Arm, Component, ESS_butterfly, Guide_gravity, Progress_bar


@dataclass(kw_only=True, repr=False)
class Slit(Component):
    xwidth: float
    yheight: float


def test_the_component_name_is_the_class_name():
    assert Slit(xwidth=0.1, yheight=0.2).component_name == "Slit"


def test_the_name_is_unset_unless_given():
    assert Slit(xwidth=0.1, yheight=0.2).name is None
    assert Slit(name="slit", xwidth=0.1, yheight=0.2).name == "slit"


def test_the_placement_defaults_to_the_absolute_origin():
    slit = Slit(xwidth=0.1, yheight=0.2)
    assert (slit.at, slit.rotated, slit.relative) == ([0, 0, 0], [0, 0, 0], None)


def test_the_mcstas_parameters_leave_out_the_name_and_the_placement():
    slit = Slit(
        name="slit",
        at=[0, 0, 1],
        rotated=[0, "omega", 0],
        relative="origin",
        xwidth=0.1,
        yheight=0.2,
    )
    assert slit.mcstas_parameters() == {"xwidth": 0.1, "yheight": 0.2}


def test_the_mcstas_parameters_leave_out_the_instance_keywords():
    slit = Slit(
        when="flag == 1",
        group="slits",
        jump="myself ITERATE 5",
        split="n_split",
        extend="hit = SCATTERED;",
        xwidth=0.1,
        yheight=0.2,
    )
    assert slit.mcstas_parameters() == {"xwidth": 0.1, "yheight": 0.2}


def test_the_instance_keywords_are_unset_by_default():
    slit = Slit(xwidth=0.1, yheight=0.2)
    assert [slit.when, slit.group, slit.jump, slit.split, slit.extend] == [None] * 5


def test_the_name_is_not_compared():
    assert Slit(name="a", xwidth=0.1, yheight=0.2) == Slit(name="b", xwidth=0.1, yheight=0.2)


def test_repr_is_a_table_of_the_fields_under_the_name():
    text = repr(Slit(name="slit", at=[0, 0, 1], xwidth=0.1, yheight=1 / 3))
    assert text.startswith("slit (Slit):\n")
    assert "rotated" in text and "[0, 0, 0]" in text
    assert "0.333333" in text
    assert repr(Slit(xwidth=0.1, yheight=0.2)).startswith("Slit:\n")


def test_repr_shows_the_other_instance_keywords_only_when_set():
    assert "split" not in repr(Slit(xwidth=0.1, yheight=0.2))
    assert "split" in repr(Slit(split=10, xwidth=0.1, yheight=0.2))
    assert "relative" in repr(Slit(xwidth=0.1, yheight=0.2))


def test_progress_bar():
    assert Progress_bar().percent == 10  # McStas's default
    origin = Progress_bar(name="origin", percent=10)
    assert origin.component_name == "Progress_bar"
    assert origin.name == "origin"
    assert origin.at == [0, 0, 0]
    assert origin.mcstas_parameters() == {"percent": 10}


def ess_butterfly(**values):
    return ESS_butterfly(
        **{
            "name": "source",
            "sector": "W",
            "beamline": 7,
            "acc_power": 2,
            "yheight": 0.03,
            "cold_frac": 0.3,
            "dist": 2,
            "focus_xw": 0.095,
            "focus_yh": 0.035,
            "n_pulses": 1,
            "Lmin": 0.1,
            "Lmax": "lambda_max",
        }
        | values
    )


def test_ess_butterfly():
    assert ess_butterfly().mcstas_parameters() == {
        "Lmin": 0.1,
        "Lmax": "lambda_max",
        "sector": "W",
        "beamline": 7,
        "acc_power": 2,
        "yheight": 0.03,
        "cold_frac": 0.3,
        "dist": 2,
        "focus_xw": 0.095,
        "focus_yh": 0.035,
        "n_pulses": 1,
        "c_performance": 1,
        "t_performance": 1,
        "target_index": 0,
        "tmax_multiplier": 3,
        "tfocus_dist": 0,
        "tfocus_time": 0,
        "tfocus_width": 0,
    }


def test_ess_butterfly_needs_only_the_wavelength_band():
    source = ESS_butterfly(Lmin=0.1, Lmax=10.0)
    assert (source.sector, source.beamline, source.acc_power) == ("W", 7, 2)
    assert source.cold_frac == 0.5  # McStas's default
    with pytest.raises(TypeError, match="Lmax"):
        ESS_butterfly(Lmin=0.1)


@pytest.mark.parametrize("sector", ['"W"', "w", "X"])
def test_ess_butterfly_needs_a_sector_letter(sector):
    source = ess_butterfly(sector=sector)
    with pytest.raises(ValueError, match="is not 'N', 'S', 'E' or 'W'"):
        source.mcstas_parameters()


def test_arm_has_no_mcstas_parameters():
    arm = Arm(name="sample_arm", at=[0, 0, 160], rotated=[0, "a3", 0], relative="origin")
    assert arm.component_name == "Arm"
    assert arm.mcstas_parameters() == {}
    assert repr(arm).startswith("sample_arm (Arm):\n")


def guide_gravity(**values):
    return Guide_gravity(
        name="guide_1",
        w1=0.089955,
        h1=0.029992,
        w2=0.082455,
        h2=0.034937,
        l=0.9994,
        mleft=1.5,
        mright=1.5,
        mtop=4,
        mbottom=4,
        **values,
    )


def test_guide_gravity():
    parameters = guide_gravity().mcstas_parameters()
    assert len(parameters) == 23
    assert parameters["l"] == 0.9994
    assert parameters["mtop"] == 4
    assert parameters["m"] == 1.0  # McStas's default


def test_guide_gravity_needs_only_its_entrance_and_length():
    guide = Guide_gravity(w1=0.05, h1=0.1, l=2)
    assert (guide.w2, guide.h2, guide.mleft) == (0, 0, -1)  # McStas's defaults
    with pytest.raises(TypeError, match="'l'"):
        Guide_gravity(w1=0.05, h1=0.1)
