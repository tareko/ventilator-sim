"""The bedside shell: parsing, commands, and the non-interactive entry point."""

import json

import pytest

from ventsim.cli import (Session, do_command, grid_from_tokens, parse_assignments,
                         parse_duration, main)


def test_durations():
    assert parse_duration("90s") == pytest.approx(1.5)
    assert parse_duration("10m") == pytest.approx(10.0)
    assert parse_duration("10") == pytest.approx(10.0)
    assert parse_duration("2h") == pytest.approx(120.0)
    assert parse_duration(" 1.5h ") == pytest.approx(90.0)
    for bad in ("", "later", "10 x", "abc"):
        with pytest.raises(ValueError):
            parse_duration(bad)


def test_assignments_land_in_the_right_namespace():
    (a,) = parse_assignments(["peep=12"])
    assert a.key == "vent.peep" and a.value == 12.0
    (a,) = parse_assignments(["vt=6"])
    assert a.key == "vent.vt_ml_kg"
    (a,) = parse_assignments(["hgb=11"])
    assert a.key == "patient.hgb"
    (a,) = parse_assignments(["mode=psv"])
    assert a.value == "psv"
    (a,) = parse_assignments(["ti=0"])
    assert a.value is None                       # 0 means "let I:E decide"
    (a,) = parse_assignments(["n_units=24"])
    assert a.value == 24 and isinstance(a.value, int)


def test_assignment_errors_are_useful():
    with pytest.raises(ValueError):
        parse_assignments(["peep"])
    with pytest.raises(ValueError):
        parse_assignments(["nonsense=1"])
    with pytest.raises(ValueError):
        parse_assignments(["mode=x"])
    with pytest.raises(ValueError):
        parse_assignments(["peep=5"], patient_only=True)     # ventilator knob
    parse_assignments(["hco3=28"], patient_only=True)        # fine


def test_search_grid_aliases():
    assert grid_from_tokens(["vt=4,6"]) == {"vt_ml_kg": [4.0, 6.0]}
    assert grid_from_tokens(["peep=5,10", "fio2=0.5"]) == {
        "peep": [5.0, 10.0], "fio2": [0.5]}
    with pytest.raises(ValueError):
        grid_from_tokens(["compliance=40"])


@pytest.fixture()
def session() -> Session:
    s = Session()
    s.start("ards_mod")
    return s


def test_session_start_and_reset(session: Session):
    assert session.scenario.key == "ards_mod"
    assert session.sim.vent.peep == 10
    session.apply(parse_assignments(["peep=18"]))
    assert session.sim.vent.peep == 18
    session.reset()
    assert session.sim.vent.peep == 10
    assert session.sim.t_min == pytest.approx(0.0)


def test_session_moves_the_blood_gas(session: Session):
    before = session.sim.current().paco2
    session.apply(parse_assignments(["rr=28"]))
    session.sim.run(20, sample_every=10, dt=0.5)
    assert session.sim.current().paco2 < before - 3


def test_panel_and_flags(session: Session):
    text = session.abg_text()
    assert "PaCO2" in text and "Pplat" in text
    assert "flags" in text
    assert "heading" in text


def test_commands(session: Session):
    out = do_command(session, "peep=14 fio2=0.7")
    assert "set" in out and "PaCO2" in out
    assert session.sim.vent.peep == 14 and session.sim.vent.fio2 == 0.7
    out = do_command(session, "vq")
    assert out.count("\n") >= 24 and "perfusion by V/Q bin" in out
    detail = do_command(session, "vq 3")
    assert "unit 3" in detail and "V/Q" in detail
    assert do_command(session, "vq 99").startswith("  there is no unit 99")
    assert len(do_command(session, "scenarios").strip().split("\n")) == 12
    assert "unknown command" in do_command(session, "frobnicate")
    assert "healthy" in do_command(session, "scenario normal").lower()


def test_clock_commands(session: Session):
    do_command(session, "run 10m")
    assert session.sim.t_min == pytest.approx(10.0, abs=1e-6)
    table = do_command(session, "watch 20m every 5m")
    rows = [ln for ln in table.split("\n") if ln.strip() and not ln.strip().startswith("t")]
    assert len(rows) == 4
    assert len(session.history) >= 5
    trend = do_command(session, "trend 3")
    assert len(trend.strip().split("\n")) == 4        # header plus three rows


def test_optimize_command_runs(session: Session):
    out = do_command(session, "optimize peep=5,10,15 fio2=0.5,0.8")
    assert "ranked" in out and "cost" in out
    assert "current" in out


def test_list_scenarios(capsys):
    assert main(["--list-scenarios"]) == 0
    printed = capsys.readouterr().out
    assert "ards_mod" in printed and "weaning" in printed


def test_json_output(capsys):
    assert main(["--scenario", "copd", "--set", "rr=20", "--run", "10m",
                 "--json", "--quiet"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["scenario"] == "copd"
    assert payload["vent"]["rr"] == 20
    assert payload["t_min"] == pytest.approx(10.0, abs=1e-6)
    assert len(payload["units"]) == 24
    assert payload["history"][-1]["paco2"] < 45


def test_bad_scenario_exits_nonzero(capsys):
    assert main(["--scenario", "not_a_case"]) == 2
    assert "unknown scenario" in capsys.readouterr().err
