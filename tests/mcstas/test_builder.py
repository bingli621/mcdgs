import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pytest

from mcdgs.core.base import Instrument
from mcdgs.core.components import Arm, ESS_butterfly, Monitor_nD, Progress_bar
from mcdgs.core.params import InstrParam
from mcdgs.instruments.trex.trex import TRex
from mcdgs.mcstas.builder import build_mcstas, default_input_dir
from mcdgs.samples import aluminum, la2ni7


@dataclass(frozen=True)
class _NoParams:
    pass


class _Beamline(Instrument):
    """Another instrument: the builder takes it unchanged."""

    params_class = _NoParams
    name = "beamline"

    def __init__(self, params=None):
        super().__init__(params)
        self.origin = Progress_bar(at=[0, 0, 0])
        self.axis = Arm(at=[0, 0, 1], rotated=[0, "2.0*omega", 0], relative="origin")
        self.monitor = Monitor_nD(xwidth=0.1, yheight=0.1, at=[0, 0, 0.5], relative="axis")

    @property
    def instr_params(self):
        return [InstrParam("omega", 30.0, "deg", "sample rotation")]


def test_components_are_added_in_order_with_their_mcstas_names(tmp_path):
    trex = TRex(sample=aluminum())
    inst = build_mcstas(trex, input_path=tmp_path)
    assert [c.name for c in inst.component_list] == list(trex.components)
    assert [c.component_name for c in inst.component_list] == [
        "Progress_bar", "ESS_butterfly", "Arm", "Guide_gravity", "Guide_gravity", "Guide_gravity",
        "Arm", "Guide_gravity", "Pol_bender", "Guide_gravity", "Arm", "Single_crystal",
    ]  # fmt: skip


def test_the_input_folder_is_created(tmp_path):
    folder = tmp_path / "new" / "in"
    inst = build_mcstas(_Beamline(), input_path=folder)
    assert folder.is_dir() and Path(inst.input_path) == folder


def test_by_default_the_input_folder_is_per_instrument_in_the_current_folder(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    beamline = _Beamline()
    inst = build_mcstas(beamline)
    assert Path(inst.input_path) == tmp_path / "mcstas_input" / "beamline"
    assert default_input_dir(beamline, "elsewhere") == Path("elsewhere/mcstas_input/beamline")


def test_the_instrument_parameters_follow_the_params(tmp_path):
    trex = TRex()
    trex.update_params(bender_in=True)
    params = {p.name: p.value for p in build_mcstas(trex, input_path=tmp_path).parameters}
    assert params == {"bender_in": 1, "b_rot": trex.params.b_rot}


def test_a_string_instrument_parameter_is_declared_as_a_string(tmp_path, monkeypatch):
    monkeypatch.setattr(
        _Beamline, "instr_params", [InstrParam("mode", "HF", "", "resolution setting")]
    )
    inst = build_mcstas(_Beamline(), input_path=tmp_path)
    assert (inst.parameters["mode"].type, inst.parameters["mode"].value) == ("string", '"HF"')


def test_the_placement_is_as_given(tmp_path):
    beamline = _Beamline()
    inst = build_mcstas(beamline, input_path=tmp_path)
    origin, axis = inst.get_component("origin"), inst.get_component("axis")
    assert origin.AT_reference is None and origin.ROTATED_reference is None  # ABSOLUTE
    assert axis.AT_data == [0, 0, 1] and axis.AT_reference == "origin"
    assert axis.ROTATED_data == [0, "2.0*omega", 0] and axis.ROTATED_reference == "origin"
    assert beamline.axis.rotated == [0, "2.0*omega", 0]  # the component itself is not changed


def test_the_instance_keywords_are_passed_on(tmp_path):
    beamline = _Beamline()
    beamline.axis.when = "omega>0"
    beamline.axis.group = "axes"
    beamline.axis.jump = "axis ITERATE 5"
    beamline.axis.split = "omega*2"
    beamline.axis.extend = "double seen = 1;"
    axis = build_mcstas(beamline, input_path=tmp_path).get_component("axis")
    assert axis.WHEN == "WHEN (omega>0)"
    assert axis.GROUP == "axes"
    assert axis.JUMP == "axis ITERATE 5"
    assert axis.SPLIT == "omega*2"
    assert "double seen = 1;" in axis.EXTEND


def test_the_instance_keywords_of_trex_choose_guide_4(tmp_path):
    inst = build_mcstas(TRex(), input_path=tmp_path)
    assert inst.get_component("guide_4a").WHEN == "WHEN (bender_in==0)"
    assert inst.get_component("bender").WHEN == "WHEN (bender_in==1)"
    assert inst.get_component("bender").ROTATED_data == [0, "b_rot", 0]
    assert inst.get_component("guide_1").WHEN == ""  # no WHEN: it always acts


def test_a_string_parameter_is_quoted(tmp_path):
    source = build_mcstas(TRex(), input_path=tmp_path).get_component("source")
    assert source.sector == '"W"'
    assert source.beamline == 7  # a number is not quoted


def test_an_unset_string_parameter_is_mcstas_null(tmp_path):
    monitor = build_mcstas(_Beamline(), input_path=tmp_path).get_component("monitor")
    assert monitor.options == '"NULL"' and monitor.filename == '"NULL"'


def test_an_expression_is_not_quoted(tmp_path):
    beamline = _Beamline()
    beamline.monitor.xwidth = "2.0*omega"  # McStasScript takes an expression
    monitor = build_mcstas(beamline, input_path=tmp_path).get_component("monitor")
    assert monitor.xwidth == "2.0*omega"


def test_a_vector_parameter_is_a_c_initializer(tmp_path):
    bender = build_mcstas(TRex(), input_path=tmp_path).get_component("bender")
    assert bender.rTopUpPar == "{0.99,0.0219,3.02,4.0,0.003}"


def test_a_cif_file_reaches_mcstas_as_its_full_path(tmp_path):
    crystal = build_mcstas(TRex(sample=la2ni7()), input_path=tmp_path).get_component("sample")
    assert crystal.reflections.endswith('La2Ni7.cif"')
    assert Path(crystal.reflections.strip('"')).is_absolute()


def test_a_field_mcstas_does_not_know_is_left_out(tmp_path):
    crystal = build_mcstas(TRex(sample=aluminum()), input_path=tmp_path).get_component("sample")
    assert "lattice_params" not in crystal.parameter_names  # the crystal computes from it
    assert crystal.recip_cell == 1 and crystal.cz == pytest.approx(1.5516, abs=1e-4)


def test_a_computed_value_mcstas_does_not_know_raises(tmp_path, monkeypatch):
    monkeypatch.setattr(Monitor_nD, "mcstas_parameters", lambda _: {"bogus": 1})
    with pytest.raises(
        ValueError, match=r"monitor \(Monitor_nD\): no McStas parameter \['bogus'\]"
    ):
        build_mcstas(_Beamline(), input_path=tmp_path)


def test_a_wrong_sector_stops_the_build(tmp_path):
    trex = TRex()
    trex.source = ESS_butterfly(Lmin=0.1, Lmax=10, sector="Q")
    with pytest.raises(ValueError, match="source: sector 'Q'"):
        build_mcstas(trex, input_path=tmp_path)


@pytest.mark.skipif(shutil.which("mcrun") is None, reason="McStas is not on the PATH")
def test_mcstas_compiles_the_written_instrument(tmp_path):
    inst = build_mcstas(TRex(sample=aluminum()), input_path=tmp_path)
    inst.write_full_instrument()
    result = subprocess.run(  # compiled, not only translated: McStas lets a bad vector through
        ["mcrun", "-c", "-n", "0", "--no-output-files", "TRex.instr"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
