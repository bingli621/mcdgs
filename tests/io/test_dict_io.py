import pytest

from mcdgs.instruments import TRex
from mcdgs.io.dict_io import from_dict, from_overrides, merge, to_dict
from mcdgs.samples import la2ni7

EXTRA_SLIT = {
    "component": "Slit",
    "xwidth": 0.01,
    "yheight": 0.01,
    "at": [0, 0, 0.5],
    "relative": "source",
}


def _with_section_after(data, after, name, section):
    """`data` with the section `name` inserted right after the key `after`."""
    items = list(data.items())
    items.insert(list(data).index(after) + 1, (name, section))
    return dict(items)


def test_a_complete_dict_gives_the_instrument_back():
    trex = TRex(sample=la2ni7())
    trex.update_params(lambda0=4.0, bender_in=True, mode="RR")
    data = to_dict(trex)
    assert from_dict(data) == trex
    assert data == to_dict(trex)  # not changed


def test_the_dict_names_the_class_then_holds_the_params_and_the_instr_params():
    trex = TRex()
    data = to_dict(trex)
    assert list(data)[:5] == ["instrument", "name", "params", "instr_params", "origin"]
    assert (data["instrument"], data["name"]) == ("TRex", "TRex")
    assert data["params"]["lambda0"] == {"value": 2.0, "unit": "Å", "comment": "central wavelength"}
    assert {name: entry["value"] for name, entry in data["instr_params"].items()} == {
        p.name: p.value for p in trex.instr_params
    }


def test_a_section_is_the_component_then_where_it_sits():
    trex = TRex()
    section = to_dict(trex)["bender"]
    assert next(iter(section)) == "component" and section["component"] == "Pol_bender"
    assert list(section)[-4:] == ["at", "rotated", "relative", "position"]
    assert section["when"] == "bender_in==1"
    assert "name" not in section and "group" not in section  # the key; unset
    assert section["position"] == trex.positions()["bender"]


def test_what_follows_from_the_rest_is_not_read():
    data = to_dict(TRex())
    data["instr_params"]["nu_bw1"]["value"] = 99.0
    data["guide_1"]["position"] = [0, 0, 99.0]  # far downstream: no part in the order
    assert from_dict(data) == TRex()


def test_a_new_component_is_built_where_its_section_is():
    data = _with_section_after(to_dict(TRex()), "source", "extra_slit", EXTRA_SLIT)
    assert list(from_dict(data).components)[:4] == [
        "origin",
        "source",
        "extra_slit",
        "source_monitor",
    ]


def test_a_component_of_the_class_without_a_section_is_left_out():
    data = to_dict(TRex())
    del data["slit_mc"]  # no other component is placed relative to it
    trex = from_dict(data)
    assert "slit_mc" not in trex.components and trex.slit_mc is None


@pytest.mark.parametrize(
    ("relative", "message"),
    [("sorce", "extra_slit: relative 'sorce' is no component"), ("guide_1", "comes after it")],
)
def test_a_relative_must_name_an_earlier_component(relative, message):
    section = EXTRA_SLIT | {"relative": relative}
    data = _with_section_after(to_dict(TRex()), "source", "extra_slit", section)
    with pytest.raises(ValueError, match=message):
        from_dict(data)


def test_a_wrong_section_says_which_one():
    data = to_dict(TRex())
    data["slit_mc"]["component"] = "Slt"
    with pytest.raises(ValueError, match="slit_mc: no such component 'Slt'"):
        from_dict(data)
    data = to_dict(TRex())
    data["slit_mc"]["xwidht"] = data["slit_mc"].pop("xwidth")
    with pytest.raises(ValueError, match="slit_mc: .*xwidht"):
        from_dict(data)
    data = to_dict(TRex())
    data["slit_mc"] = 0.02
    with pytest.raises(ValueError, match="slit_mc: a component section is a mapping"):
        from_dict(data)


@pytest.mark.parametrize(
    ("params", "match"),
    [
        ({"lambda0": {"value": 0.4, "unit": "nm"}}, "params.lambda0: unit 'nm', but it is 'Å'"),
        ({"lamda0": 4.0}, r"params: no such param \['lamda0'\]"),
    ],
)
def test_a_wrong_param_is_refused(params, match):
    with pytest.raises(ValueError, match=match):
        from_overrides({"instrument": "TRex", "params": params})


def test_a_param_may_be_a_bare_value_and_its_comment_is_not_read():
    trex = from_overrides(
        {"instrument": "TRex", "params": {"lambda0": 4.0, "rrm": {"value": 6, "comment": "x"}}}
    )
    assert (trex.params.lambda0, trex.params.rrm) == (4.0, 6)


def test_overrides_are_laid_over_the_default_instrument():
    trex = from_overrides({"instrument": "TRex", "slit_mc": {"xwidth": 0.02, "yheight": None}})
    assert (trex.slit_mc.xwidth, trex.slit_mc.yheight) == (0.02, TRex().slit_mc.yheight)
    assert trex.guide_1 == TRex().guide_1  # left out: as made


def test_the_class_is_given_what_its_default_instrument_needs():
    sample = la2ni7()
    assert from_overrides({"instrument": "TRex"}, sample=sample) == TRex(sample=sample)


def test_a_dict_must_name_its_instrument():
    data = to_dict(TRex())
    del data["instrument"]
    with pytest.raises(ValueError, match="instrument: None is no instrument"):
        from_dict(data)
    with pytest.raises(ValueError, match="instrument: 'Nothing' is no instrument"):
        from_overrides({"instrument": "Nothing"})


@pytest.mark.parametrize("name", ["positions", "update_params", "add_components", "_params"])
def test_a_section_named_like_an_attribute_of_the_instrument_is_refused(name):
    data = to_dict(TRex())
    data[name] = {"component": "Arm", "at": [0, 0, 0]}
    with pytest.raises(ValueError, match=f"{name}: reserved name"):
        from_dict(data)


def test_a_component_named_like_a_key_of_the_file_is_refused():
    trex = TRex()
    trex.instrument = trex.origin
    with pytest.raises(ValueError, match="instrument: reserved name in a file"):
        to_dict(trex)


def test_merge_lays_one_dict_over_another():
    base = {"a": {"x": 1, "y": 2}, "b": [1, 2]}
    assert merge(base, {"a": {"x": 5, "y": None}, "b": [3]}) == {"a": {"x": 5, "y": 2}, "b": [3]}
    assert base == {"a": {"x": 1, "y": 2}, "b": [1, 2]}  # not changed


def test_merge_puts_a_new_key_after_the_key_before_it_in_the_overrides():
    defaults = {"a": 1, "b": 2, "c": 3}
    assert list(merge(defaults, {"b": 5, "new": 0})) == ["a", "b", "new", "c"]
    assert list(merge(defaults, {"new": 0, "c": 5})) == ["new", "a", "b", "c"]  # first: first
    assert list(merge(defaults, {"c": 5, "a": 5})) == ["a", "b", "c"]  # the defaults' order kept
