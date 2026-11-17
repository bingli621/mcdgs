import pytest

from mcdgs.instruments.trex.model import TRex
from mcdgs.mcstas import RunType, Simulation, Workspace


def test_a_scan_has_a_point_per_value_of_its_variable():
    points = Simulation(TRex()).get_scan_points("b_rot", range(3), {"bender_in": 1})
    assert points == [
        {"bender_in": 1, "b_rot": 0.0},
        {"bender_in": 1, "b_rot": 1.0},
        {"bender_in": 1, "b_rot": 2.0},
    ]


@pytest.mark.parametrize(("variable", "instr_params"), [("b_rto", {}), ("b_rot", {"bender": 1})])
def test_a_scan_of_an_unknown_parameter_raises(variable, instr_params):
    with pytest.raises(KeyError):
        Simulation(TRex()).get_scan_points(variable, [1, 2], instr_params)


@pytest.mark.parametrize(
    ("method", "kwargs"),
    [
        ("run_to_sample", {}),
        ("run_from_sample", {"beam_run": "beam"}),
        ("scan_from_sample", {"variable": "b_rot", "values": [1], "beam_run": "beam"}),
    ],
)
def test_a_split_run_without_sample_pos_raises_before_anything(tmp_path, method, kwargs):
    trex = TRex()
    del trex.sample_pos
    sim = Simulation(trex, Workspace(tmp_path))
    with pytest.raises(ValueError, match="sample_pos"):
        getattr(sim, method)(**kwargs)
    assert not any(tmp_path.iterdir())  # no run folder


def test_a_run_from_the_sample_without_a_run_to_the_sample_raises_before_anything(tmp_path):
    sim = Simulation(TRex(), Workspace(tmp_path))
    with pytest.raises(FileNotFoundError, match="run_to_sample first"):
        sim.run_from_sample()
    with pytest.raises(FileNotFoundError, match="no finished run 'missing'"):
        sim.scan_from_sample("b_rot", [1], beam_run="missing")
    assert not any(tmp_path.iterdir())  # no run folder


def test_the_beam_of_an_earlier_session_is_named_not_guessed(tmp_path, finished_run):
    beam = finished_run(tmp_path / "runs" / "beam", 1, RunType.RUN_TO_SAMPLE)
    finished_run(tmp_path / "runs" / "newer_beam", 2, RunType.RUN_TO_SAMPLE, instrument="other")
    sim = Simulation(TRex(), Workspace(tmp_path))
    assert sim.beam_run is None
    with pytest.raises(FileNotFoundError):
        sim.find_beam_run()  # none given, none of this simulation
    assert sim.find_beam_run("beam") == beam
    sim.beam_run = "beam"  # as run_to_sample leaves it
    assert sim.find_beam_run() == beam
    assert sim.find_beam_run("newer_beam").name == "newer_beam"  # the one given comes first


def test_only_a_run_to_the_sample_is_a_beam_run_and_it_is_not_replaced(tmp_path, finished_run):
    finished_run(tmp_path / "runs" / "beam", 1, RunType.RUN_TO_SAMPLE)
    finished_run(tmp_path / "runs" / "scan", 2, RunType.SCAN)
    sim = Simulation(TRex(), Workspace(tmp_path))
    with pytest.raises(ValueError, match="'scan' is a scan run, not a run_to_sample"):
        sim.run_from_sample(beam_run="scan")
    with pytest.raises(ValueError, match="folder_name 'beam'"):
        sim.run_from_sample("beam", beam_run="beam")
    assert (tmp_path / "runs" / "beam" / "run.json").is_file()  # still there


def test_the_folders_are_those_of_the_workspace(tmp_path, monkeypatch):
    monkeypatch.delenv("MCDGS_WORKSPACE", raising=False)
    monkeypatch.chdir(tmp_path)
    sim = Simulation(TRex())
    assert sim.workspace == Workspace(tmp_path / "mcdgs_workspace")
    monkeypatch.chdir(tmp_path.parent)  # fixed when the simulation was created
    assert sim.workspace.root == tmp_path / "mcdgs_workspace"
    assert Simulation(TRex(), Workspace(tmp_path / "w")).workspace.root == tmp_path / "w"


def test_the_record_of_a_run_has_every_instrument_parameter_it_used(tmp_path):
    sim = Simulation(TRex(), Workspace(tmp_path))
    run, _ = sim.prepare("mine", RunType.RUN, {"bender_in": 1, "omega": 30}, ncount=1e4, seed=3)
    defaults = {p.name: p.value for p in TRex().instr_params}
    assert run.instr_params.keys() == defaults.keys()
    assert (run.instr_params["bender_in"], run.instr_params["omega"]) == (1.0, 30.0)  # as given
    assert run.instr_params["b_rot"] == defaults["b_rot"]  # not given: as the instrument has it
    assert (run.instrument, run.ncount, run.seed) == ("TRex", 1e4, 3)
    assert run.folder == tmp_path / "runs" / "mine" and not run.folder.exists()  # not saved yet


def test_the_record_of_a_scan_has_its_variable_apart_from_the_other_parameters(tmp_path):
    sim = Simulation(TRex(), Workspace(tmp_path))
    run, _ = sim.prepare(
        None, RunType.SCAN_FROM_SAMPLE, {"bender_in": 1}, variable="omega", values=[0.0, 1.0]
    )
    assert run.folder.name.endswith("_scan_omega_from_sample")
    assert (run.variable, run.values) == ("omega", [0.0, 1.0])
    assert "omega" not in run.instr_params and run.instr_params["bender_in"] == 1.0
    assert run.started.tzinfo is not None


def test_the_record_of_a_run_from_the_sample_names_its_beam_run(tmp_path):
    sim = Simulation(TRex(), Workspace(tmp_path))
    run, _ = sim.prepare("a", RunType.RUN_FROM_SAMPLE, {}, beam_run="to_sample")
    assert run.beam_run == "to_sample" and run.ncount is None


def test_a_run_mcstasscript_starts_itself_is_prepared_without_nexus(tmp_path):
    sim = Simulation(TRex(), Workspace(tmp_path))
    run, inst = sim.prepare("diagnose", ncount=1e4, seed=3, nexus=False, mpi="auto")
    assert run.run_type is RunType.RUN and run.folder == tmp_path / "runs" / "diagnose"
    assert run.instr_params == {p.name: float(p.value) for p in TRex().instr_params}
    assert inst._run_settings["NeXus"] is False  # McStasScript keeps them private
    assert inst._run_settings["output_path"] == str(run.folder)
    assert (inst._run_settings["ncount"], inst._run_settings["seed"]) == (1e4, 3)
