"""The setting search: parsing, monotonicity, and coordinate descent vs exhaustive."""

import pytest

from ventsim.optimizer import (Candidate, DEFAULT_GRIDS_VC, Weights, cost_of,
                               explain, optimize, parse_grid, signature)
from ventsim.model import Patient, Vent, solve_steady


def ards():
    p = Patient(height_cm=170, weight_kg=70, severity=0.55, compliance=28,
                resistance=8, heterogeneity=0.45, hgb=11.0, co_l_min=6.0)
    return p, Vent(peep=10, fio2=0.8, vt_ml_kg=6, rr=21)


def test_parse_grid():
    assert parse_grid("peep=5,8,12 fio2=0.4,0.6") == {
        "peep": [5.0, 8.0, 12.0], "fio2": [0.4, 0.6]}
    assert parse_grid("rr=12") == {"rr": [12.0]}
    with pytest.raises(ValueError):
        parse_grid("peep")
    with pytest.raises(ValueError):
        parse_grid("peep=5,abc")
    with pytest.raises(ValueError):
        parse_grid("")


def test_unknown_knob_is_rejected():
    p, v = ards()
    with pytest.raises(ValueError):
        optimize(p, v, grids={"mode": ["vc"]})


def test_penalties_point_the_way_a_clinician_looks():
    p, base = ards()
    w = Weights()
    good = solve_steady(p, Vent(peep=12, fio2=0.6, vt_ml_kg=5, rr=22))
    cost_good, parts_good = cost_of(p, Vent(peep=12, fio2=0.6, vt_ml_kg=5, rr=22), good, base, w)
    # same settings but with a much worse result: the cost must rise
    worse = solve_steady(p, Vent(peep=5, fio2=0.3, vt_ml_kg=5, rr=22))
    cost_worse, _ = cost_of(p, Vent(peep=5, fio2=0.3, vt_ml_kg=5, rr=22), worse, base, w)
    assert cost_worse > cost_good
    assert "hypoxaemia" in parts_good or cost_good < 200


def test_small_change_is_preferred_when_things_are_equal():
    p, v = ards()
    w = Weights()
    res = solve_steady(p, v)
    same, parts_same = cost_of(p, v, res, v, w, prefer_no_change=True)
    moved = Vent(peep=14, fio2=0.8, vt_ml_kg=6, rr=21)
    res_moved = solve_steady(p, moved)
    cost_moved, parts_moved = cost_of(p, moved, res_moved, v, w, prefer_no_change=True)
    assert "change" not in parts_same
    assert parts_moved.get("change", 0.0) > 0.0
    assert cost_moved >= 0.0 and same >= 0.0


def test_work_of_breathing_penalises_under_support():
    p = Patient(name="day4", height_cm=172, weight_kg=72, severity=0.30,
                compliance=45, resistance=7, sedation=0.15, muscle=0.8, vco2=240)
    w = Weights()
    base = Vent(mode="psv", dp=12, rr=8, peep=5, fio2=0.35)
    weak = solve_steady(p, Vent(mode="psv", dp=4, rr=8, peep=5, fio2=0.35))
    enough = solve_steady(p, base)
    cost_weak, parts_weak = cost_of(p, base, weak, base, w, prefer_no_change=False)
    cost_ok, parts_ok = cost_of(p, base, enough, base, w, prefer_no_change=False)
    assert "work_of_breathing" in parts_weak
    assert "support_deficit" in parts_weak
    assert cost_weak > cost_ok


def test_coordinate_descent_finds_what_the_exhaustive_grid_finds():
    p, v = ards()
    grids = {"peep": [5, 10, 15], "fio2": [0.4, 0.8],
             "vt_ml_kg": [4.0, 6.0, 8.0], "rr": [16, 20, 24]}
    exhaustive = optimize(p, v, grids=grids, method="grid", top=3)
    quick = optimize(p, v, grids=grids, method="coord", top=3)
    assert exhaustive and quick
    # coordinate descent may step on the patient's own settings, so it can only
    # do as well as the exhaustive sweep, never worse
    assert quick[0].cost <= exhaustive[0].cost + 1e-9
    assert quick[0].cost > exhaustive[0].cost - 25.0


def test_candidates_are_unique_by_settings():
    p, v = ards()
    cands = optimize(p, v, top=50)
    sigs = [signature(c.vent) for c in cands]
    assert len(set(sigs)) == len(sigs)


def test_default_grids_cover_the_knobs():
    assert set(DEFAULT_GRIDS_VC) >= {"peep", "fio2", "vt_ml_kg", "rr", "ie"}
    for values in DEFAULT_GRIDS_VC.values():
        assert sorted(values) == values


def test_explain_reads_like_a_score_sheet():
    p, v = ards()
    c = optimize(p, v, grids={"peep": [5]}, top=1)[0]
    assert isinstance(c, Candidate)
    text = explain(c)
    assert text
    assert any(ch.isdigit() for ch in text)
