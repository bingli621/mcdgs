from dataclasses import dataclass

import pytest

from mcdgs.core.components import (
    COMPONENTS,
    Arm,
    Component,
    ESS_butterfly,
    Guide_gravity,
    Incoherent,
    Monitor_nD,
    NXdisk_chopper,
    Phonon_simple,
    Pol_bender,
    Polygon_ESS_butterfly,
    Progress_bar,
)


@dataclass(kw_only=True, repr=False)
class Slit(Component):
    xwidth: float
    yheight: float


def test_the_component_name_is_the_class_name():
    assert Slit(xwidth=0.1, yheight=0.2).component_name == "Slit"


def test_every_component_class_is_one_mcstas_component():
    import mcdgs.instruments.trex.model  # noqa: F401  (defines its components: none subclassed)

    # only the class of a McStas component that extends another subclasses that one
    extends = {Polygon_ESS_butterfly: ESS_butterfly}
    assert all(c.__bases__ == (extends.get(c, Component),) for c in COMPONENTS.values())
    assert {
        c.__name__: [s.__name__ for s in c.__subclasses__()]
        for c in COMPONENTS.values()
        if c.__subclasses__()
    } == {"ESS_butterfly": ["Polygon_ESS_butterfly"]}


def test_by_default_a_component_is_at_the_absolute_origin_with_no_instance_keyword():
    slit = Slit(xwidth=0.1, yheight=0.2)
    assert (slit.at, slit.rotated, slit.relative) == ([0, 0, 0], [0, 0, 0], None)
    assert [slit.when, slit.group, slit.jump, slit.split, slit.extend] == [None] * 5


def test_the_mcstas_parameters_leave_out_the_instance_keywords():
    slit = Slit(
        name="slit",
        at=[0, 0, 1],
        rotated=[0, "omega", 0],
        relative="origin",
        when="flag == 1",
        group="slits",
        jump="myself ITERATE 5",
        split="n_split",
        extend="hit = SCATTERED;",
        xwidth=0.1,
        yheight=0.2,
    )
    assert slit.mcstas_parameters() == {"xwidth": 0.1, "yheight": 0.2}


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
    assert Progress_bar(percent=20).mcstas_parameters() == {"percent": 20}


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


def polygon_ess_butterfly(**values):
    """A Polygon_ESS_butterfly with one disk, with ``values`` replacing its defaults."""
    disk = {
        "speed": "nu_ps1",
        "delay": "delay_ps1",
        "beam": 180,
        "edges": [-10, 10, 107.5, 142.5],
        "path": 107.95,
        "aperture": 3.2,
    }
    return Polygon_ESS_butterfly(
        **{"name": "source", "Lmin": 1, "Lmax": 3, "choppers": [disk]} | values
    )


def test_polygon_ess_butterfly_has_the_parameters_of_ess_butterfly_and_its_own():
    parameters = polygon_ess_butterfly().mcstas_parameters()
    assert ess_butterfly().mcstas_parameters().keys() < parameters.keys()
    assert parameters.keys() - ess_butterfly().mcstas_parameters().keys() == {
        "choppers",
        "chopper_count",
        "filename",
        "path_spread_fraction",
        "noise_fraction",
        "emission_depth",
        "use_region",
        "resample",
        "save_polygons",
        "verify_acceptance",
    }
    assert parameters["chopper_count"] == 1
    assert parameters["filename"] == "source"  # unset: the component's name, not McStas's NULL


def test_polygon_ess_butterfly_has_the_defaults_of_ess_butterfly():
    masked, plain = polygon_ess_butterfly(), ESS_butterfly(name="source", Lmin=1, Lmax=3)
    assert {
        k: v for k, v in masked.mcstas_parameters().items() if k in plain.mcstas_parameters()
    } == (plain.mcstas_parameters())
    assert (masked.frequency, masked.pulse_length) == (plain.frequency, plain.pulse_length)


def test_polygon_ess_butterfly_without_choppers_is_ess_butterfly():
    source = Polygon_ESS_butterfly(Lmin=1, Lmax=3)
    assert source.choppers == []
    assert source.mcstas_parameters()["chopper_count"] == 0


@pytest.mark.parametrize(
    ("values", "reason"),
    [
        ({"tfocus_width": 1e-4}, "time focusing"),
        ({"use_region": 0}, "use_region=0"),
        ({"noise_fraction": 0.1}, "noise_fraction=0.1"),
    ],
)
def test_polygon_ess_butterfly_runs_without_resample_where_mcstas_refuses_it(
    values, reason, capsys
):
    source = polygon_ess_butterfly(**values)
    parameters = source.mcstas_parameters()
    assert f"resample=1 does not go with {reason};" in capsys.readouterr().out
    assert parameters["resample"] == 0
    assert source.resample == 1  # the field stays as it is


def test_without_resample_polygon_ess_butterfly_takes_time_focusing_silently(capsys):
    source = polygon_ess_butterfly(resample=0, tfocus_width=1e-4, noise_fraction=0.1)
    assert source.mcstas_parameters()["resample"] == 0
    assert capsys.readouterr().out == ""


def test_polygon_ess_butterfly_names_every_reason_it_runs_without_resample(capsys):
    polygon_ess_butterfly(use_region=0, noise_fraction=0.1).mcstas_parameters()
    assert "with use_region=0, noise_fraction=0.1;" in capsys.readouterr().out


def test_arm_has_no_mcstas_parameters():
    arm = Arm(name="sample_arm", at=[0, 0, 160], rotated=[0, "a3", 0], relative="origin")
    assert arm.mcstas_parameters() == {}


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


def test_monitor_nd():
    monitor = Monitor_nD(
        name="lambda_monitor",
        at=[0, 0, 0.5],
        relative="sample_arm",
        xwidth=0.1,
        yheight=0.1,
        options="lambda limits=[0.5 5] bins=50",
    )
    parameters = monitor.mcstas_parameters()
    assert len(parameters) == 19
    assert parameters["options"] == "lambda limits=[0.5 5] bins=50"
    assert "user0" not in parameters and "username0" not in parameters


def test_monitor_nd_defaults_are_mcstas_s():
    monitor = Monitor_nD()
    assert (monitor.min, monitor.max) == (-1e40, 1e40)
    assert (monitor.options, monitor.filename, monitor.geometry) == (None, None, None)
    assert (monitor.bins, monitor.restore_neutron, monitor.nowritefile) == (0, 0, 0)


def test_monitor_nd_bin_limits_reads_every_axis_of_options():
    monitor = Monitor_nD(options="time limits=[0.01 0.1] bins=1800 lambda limits=[0.1 10] bins=495")
    assert monitor.bin_limits == {"time": (0.01, 0.1, 1800), "lambda": (0.1, 10.0, 495)}
    assert Monitor_nD(options="lambda limits=[0.5 5] bins=50").bin_limits == {
        "lambda": (0.5, 5.0, 50)
    }


def test_monitor_nd_bin_limits_is_empty_without_axes():
    assert Monitor_nD().bin_limits == {}
    assert Monitor_nD(options="x bins=4, y bins=3, list all, pixel id, t").bin_limits == {}


def test_a_bender_keeps_the_defaults_of_mcstas():
    bender = Pol_bender(xwidth=0.06, yheight=0.04, length=0.05, radius=7.2)
    parameters = bender.mcstas_parameters()
    assert (parameters["G"], parameters["nslit"], parameters["drawOption"]) == (9.8, 1, 1)
    assert parameters["rRightDownPar"] == [0.99, 0.0219, 6.07, 2.0, 0.003]


def test_a_bender_needs_its_size():
    with pytest.raises(TypeError):
        Pol_bender(xwidth=0.06, yheight=0.04, length=0.05)


def test_two_benders_share_no_mirror():
    first = Pol_bender(xwidth=0.06, yheight=0.04, length=0.05, radius=7.2)
    second = Pol_bender(xwidth=0.06, yheight=0.04, length=0.05, radius=7.2)
    first.rTopUpPar[3] = 4.0
    assert second.rTopUpPar[3] == 2.0
    assert first.rTopDownPar[3] == 2.0


def test_the_chopper_defaults_are_those_of_mcstas():
    parameters = NXdisk_chopper(slit_edges=[-10, 10]).mcstas_parameters()
    assert parameters == {
        "slit_edges": [-10, 10],
        "n_edges": 2,
        "radius": 0.35,
        "yheight": 0,
        "xwidth": 0,
        "nu": 0,
        "delay": 0,
        "park_angle": 0,
        "beam_angle": 0,
        "zero_angle": 0,
        "jitter": 0,
        "abs_out": 0,
        "verbose": 0,
    }


def test_the_chopper_counts_its_edges():
    chopper = NXdisk_chopper(slit_edges=[0, 10, 90, 100])
    assert chopper.mcstas_parameters()["n_edges"] == 4


def test_the_chopper_requires_its_edges():
    with pytest.raises(TypeError):
        NXdisk_chopper()


PHONON = {  # what Phonon_simple requires, as in McStas
    "radius": 0.005, "yheight": 0.01, "sigma_abs": 0.171, "sigma_inc": 0.003,
    "a": 4.9508, "b": 9.405, "M": 207.2, "c": 13.75, "DW": 1.0, "T": 300.0,
}  # fmt: skip


def test_the_phonon_sample_defaults_are_those_of_mcstas():
    assert Phonon_simple(**PHONON).mcstas_parameters() == PHONON | {
        "target_x": 0,
        "target_y": 0,
        "target_z": 0,
        "target_index": 0,
        "focus_r": 0,
        "focus_xw": 0,
        "focus_yh": 0,
        "focus_aw": 0,
        "focus_ah": 0,
        "gap": 0,
        "e_steps_low": 50,
        "e_steps_high": 50,
    }


@pytest.mark.parametrize("missing", sorted(PHONON))
def test_the_phonon_sample_requires_what_mcstas_requires(missing):
    with pytest.raises(TypeError):
        Phonon_simple(**{k: v for k, v in PHONON.items() if k != missing})


def test_the_phonon_branch_ends_at_four_times_the_sound_velocity_over_the_lattice_constant():
    assert Phonon_simple(**PHONON).branch_maximum == pytest.approx(4 * 13.75 / 4.9508)


def test_the_incoherent_sample_defaults_are_those_of_mcstas():
    assert Incoherent().mcstas_parameters() == {
        "geometry": None,
        "radius": 0,
        "xwidth": 0,
        "yheight": 0,
        "zdepth": 0,
        "thickness": 0,
        "target_x": 0,
        "target_y": 0,
        "target_z": 0,
        "target_index": 0,
        "focus_r": 0,
        "focus_xw": 0,
        "focus_yh": 0,
        "focus_aw": 0,
        "focus_ah": 0,
        "pack": 1,
        "p_interact": 1,
        "f_QE": 0,
        "gamma": 0,
        "Etrans": 0,
        "deltaE": 0,
        "sigma_abs": 5.08,
        "sigma_inc": 5.08,
        "Vc": 13.827,
        "concentric": 0,
        "order": 0,
    }


def test_the_incoherent_scattering_length_follows_the_density_of_atoms():
    assert Incoherent().scattering_length == pytest.approx(13.827 / 508)  # m
    assert Incoherent(pack=0.5).scattering_length == pytest.approx(2 * 13.827 / 508)
    assert Incoherent(sigma_inc=80.27, Vc=30.0).scattering_length == pytest.approx(30 / 8027)
