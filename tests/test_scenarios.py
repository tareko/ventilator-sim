"""The teaching cases must build, and must teach what they claim to teach."""

import pytest

from ventsim.scenarios import BY_KEY, SCENARIOS, get


def test_scenario_registry_is_consistent():
    assert len(SCENARIOS) == 12
    assert set(BY_KEY) == {s.key for s in SCENARIOS}
    assert get("normal").key == "normal"
    for s in SCENARIOS:
        assert s.title and s.story and s.expect
        assert s.try_, f"{s.key} has no suggested experiments"
        assert s.vent.mode in ("vc", "pc", "psv")
        assert 0.21 <= s.vent.fio2 <= 1.0


def test_every_scenario_builds_and_reports():
    for s in SCENARIOS:
        sim = s.build()
        r = sim.current()
        assert r.paco2 > 5 and r.paco2 < 160, s.key
        assert r.pao2 > 20, s.key
        assert 7.05 < r.ph < 7.60, s.key
        assert r.vt_ml > 0, s.key


def test_only_the_chronic_case_needs_renal_adaptation():
    chronic = [s.key for s in SCENARIOS if s.adapt_hours]
    assert chronic == ["chronic_hypercapnic"]
    assert BY_KEY["chronic_hypercapnic"].adapt_hours >= 24


def test_the_healthy_case_is_the_ruler():
    r = get("normal").build().current()
    assert r.paco2 == pytest.approx(41.5, abs=3.0)
    assert r.pf > 380
    assert r.shunt_pct < 5
    assert r.pplat < 18
    assert not [lvl for lvl, _ in r.flags if lvl == "alert"]


def test_ards_severity_orders_by_mechanics_not_by_pf():
    mild = get("ards_mild").build().current()
    mod = get("ards_mod").build().current()
    severe = get("ards_severe").build().current()
    assert mild.pplat < mod.pplat < severe.pplat
    assert mild.vd_vt < mod.vd_vt < severe.vd_vt
    assert severe.dp_stat > 15 and severe.pplat > 35
    assert mod.shunt_pct > 15


def test_pneumonia_is_a_shunt_problem_that_peep_partly_fixes():
    s = get("pneumonia")
    sim = s.build()
    r = sim.current()
    assert r.shunt_pct > 25
    assert r.pf < 120
    better = sim.set_vent(peep=18) and sim.target()
    assert better.shunt_pct < r.shunt_pct


def test_embolism_is_dead_space_not_shunt():
    r = get("pe").build().current()
    assert r.vd_vt > 0.40
    assert r.shunt_pct < 12


def test_asthma_traps_air_and_limits_flow():
    r = get("asthma").build().current()
    assert r.auto_peep > 15
    assert r.co_l_min < 4.0
    assert r.pplat > 35


def test_copd_can_ventilate_but_aims_high():
    s = get("copd")
    sim = s.build()
    baseline = sim.current()
    assert baseline.auto_peep > 0.5
    # a longer exhalation should reduce the trapping
    sim.set_vent(ie=3.0)
    slower = sim.target()
    assert slower.auto_peep < baseline.auto_peep


def test_fibrosis_breathes_fast_and_low_compliance():
    r = get("fibrosis").build().current()
    assert r.paco2 < 33
    assert r.vd_vt > 0.30
    assert r.recruit_pct > 90


def test_obesity_splits_the_pressure_between_chest_and_lung():
    r = get("obesity").build().current()
    assert r.pplat > 30
    assert r.shunt_pct > 15


def test_chronic_case_arrives_compensated():
    r = get("chronic_hypercapnic").build().current()
    assert r.paco2 > 50
    assert r.hco3 > 28
    assert r.be > 2
    assert 7.33 < r.ph < 7.40


def test_weaning_case_shows_effort():
    r = get("weaning").build().current()
    assert r.mode == "psv"
    assert r.effort_ratio < 0.95
    assert r.drive > 1.1
