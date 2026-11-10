"""Compiled checks for the train-aware Monitor_nD."""

from test_tof_train_chopper import INSTRUMENT, _run


def test_time_axis_auto_splits_train_samples(tmp_path, monkeypatch):
    import test_tof_train_chopper as chopper

    monitor = """COMPONENT detector = Monitor_nD(
    options="t limits=[0 1] bins=10", filename="detector.dat",
    xwidth=0.08, yheight=0.08, restore_neutron=1, adaptive_target=1)"""
    original = """COMPONENT detector = PSD_monitor(
    nx=1, ny=1, filename="detector.dat", xwidth=0.08, yheight=0.08,
    restore_neutron=1)"""
    instrument = INSTRUMENT.replace(original, monitor).replace("target_tsplit=0", "target_tsplit=7")
    # Exercise coexistence with a stock component that imports McStas's own
    # monitor_nd-lib, rather than only testing the private library in isolation.
    instrument = instrument.replace(
        "FINALLY\n",
        '''COMPONENT stock_monitor = Monitor_nD_noacc(
    options="x bins=10", filename="stock.dat", xwidth=0.08, yheight=0.08)
AT (0, 0, 1.1) ABSOLUTE

FINALLY
''',
    )
    monkeypatch.setattr(chopper, "INSTRUMENT", instrument)
    assert _run(tmp_path / "auto", nu=14) > 0
    monkeypatch.setattr(
        chopper,
        "INSTRUMENT",
        instrument.replace("adaptive_target=1", "adaptive_target=1, tsplit=0"),
    )
    assert _run(tmp_path / "unsplit", nu=14) > 0

    def count(folder):
        lines = (tmp_path / folder / "run/detector.dat").read_text().splitlines()
        return float(next(line for line in lines if line.startswith("# values:")).split()[4])

    assert count("auto") > 3 * count("unsplit")
    assert (tmp_path / "auto/run/stock.dat").is_file()
