"""Contract tests for the opt-in experimental TOF train integration."""

from pathlib import Path

import pytest

from mcdgs.instruments.trex.model import TRex
from mcdgs.mcstas.builder import build_mcstas
from mcdgs.mcstas.simulation import Simulation


def test_train_build_uses_local_overrides_and_persistent_feedback(tmp_path):
    inst = build_mcstas(TRex(), input_path=tmp_path, tof_train=True)
    assert inst.input_path == str(tmp_path)
    inst.write_full_instrument()
    text = (tmp_path / "TRex.instr").read_text()
    assert "tof_train_components" in text
    component_dir = Path(inst.component_reader.component_path["ESS_butterfly"]).parent
    assert (component_dir / "trex-tof-train-butterfly-geometry.c").is_file()
    assert Path(inst.component_reader.component_path["NXdisk_chopper"]).parent == component_dir
    assert not (tmp_path / "ESS_butterfly.comp").exists()
    params = {p.name: p.value for p in inst.parameters}
    assert params["N_trains_par"] == 300
    assert params["target_tsplit"] == 10
    assert "long * train_total_N_sent" in text
    assert "train_total_N_sent = &tof_train_sent" in text
    assert "train_total_arrived = &tof_train_arrived" in text
    assert "N_trains = N_trains_par" in text
    assert text.index("COMPONENT mc2") < text.index("COMPONENT pre_sample_monitor")
    assert inst.get_component("pre_sample_monitor").adaptive_target == 1
    for name in ("source_monitor", "beam_monitor", "beamstop_monitor", "sample_monitor"):
        assert inst.get_component(name).adaptive_target == 0
    assert inst.get_component("source").target_tsplit == "target_tsplit"


def test_train_warns_on_particle_splitting(tmp_path):
    trex = TRex()
    trex.mc1.split = "2"
    with pytest.warns(UserWarning, match="after the last chopper.*absorbing"):
        inst = build_mcstas(trex, input_path=tmp_path, tof_train=True)
    assert inst.get_component("mc1").SPLIT == "2"


def test_train_still_rejects_jumps(tmp_path):
    trex = TRex()
    trex.mc1.jump = "source"
    with pytest.raises(ValueError, match="jumps"):
        build_mcstas(trex, input_path=tmp_path, tof_train=True)


def test_train_does_not_shadow_or_overwrite_local_components(tmp_path):
    local = tmp_path / "ESS_butterfly.comp"
    local.write_text("user component")
    with pytest.raises(FileExistsError, match="shadowed"):
        build_mcstas(TRex(), input_path=tmp_path, tof_train=True)
    assert local.read_text() == "user component"


def test_train_requires_mcpl_output_override_before_replacing_output(tmp_path, monkeypatch):
    original_is_file = Path.is_file
    monkeypatch.setattr(
        Path, "is_file", lambda p: False if p.name == "MCPL_output.comp" else original_is_file(p)
    )
    output = tmp_path / "mcstas_output" / "previous"
    output.mkdir(parents=True)
    sentinel = output / "keep"
    sentinel.write_text("original")
    sim = Simulation(TRex(), work_dir=tmp_path, tof_train=True)
    with pytest.raises(FileNotFoundError, match="MCPL output override missing"):
        sim.run_to_sample(folder_name="previous")
    assert sentinel.read_text() == "original"


def test_train_rejects_shadowing_mcpl_output(tmp_path):
    component = (
        Path(__file__).parents[2] / "src/mcdgs/instruments/trex/"
        "tof_train_components/MCPL_output.comp"
    )
    if not component.is_file():
        pytest.skip("experimental MCPL_output component not yet available")
    sim = Simulation(TRex(), work_dir=tmp_path, tof_train=True)
    sim.input_dir.mkdir(parents=True)
    local = sim.input_dir / "MCPL_output.comp"
    local.write_text("user component")
    with pytest.raises(FileExistsError, match="shadowed"):
        sim.run_to_sample()
    assert local.read_text() == "user component"


def test_train_rejects_custom_extend_and_preserves_existing_output(tmp_path):
    trex = TRex()
    trex.mc1.extend = "t += 1;"
    output = tmp_path / "output" / "previous"
    output.mkdir(parents=True)
    sentinel = output / "keep"
    sentinel.write_text("original")
    sim = Simulation(trex, output_dir=output.parent, work_dir=tmp_path, tof_train=True)
    with pytest.raises(ValueError, match="custom EXTEND"):
        sim.run(folder_name="previous")
    assert sentinel.read_text() == "original"


def test_train_rejects_moved_control_monitor(tmp_path):
    trex = TRex()
    trex.pre_sample_monitor.at = [1, 0, 0]
    with pytest.raises(ValueError, match="monitor placement"):
        build_mcstas(trex, input_path=tmp_path, tof_train=True)
