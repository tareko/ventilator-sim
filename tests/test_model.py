"""Steady-state behaviour: mechanics, gas exchange, flags, spontaneous breathing."""

import pytest

from ventsim import Lung, Patient, Vent, solve_steady
from ventsim.model import result_flags


def pats(**over) -> Patient:
    base = dict(name="adult", height_cm=172, weight_kg=75, compliance=60, resistance=5)
    base.update(over)
    return Patient(**base)


def ards(**over) -> Patient:
    base = dict(height_cm=170, weight_kg=70, severity=0.55, compliance=28,
                resistance=8, heterogeneity=0.45, hgb=11.0, co_l_min=6.0)
    base.update(over)
    return Patient(**base)


def test_ventilation_clears_co2_in_the_expected_direction():
    p = pats()
    low = solve_steady(p, Vent(rr=8, vt_ml_kg=8))
    mid = solve_steady(p, Vent(rr=12, vt_ml_kg=8))
    high = solve_steady(p, Vent(rr=24, vt_ml_kg=8))
    assert low.paco2 > mid.paco2 > high.paco2
    # roughly halving PaCO2 for roughly doubling alveolar ventilation
    assert high.paco2 == pytest.approx(mid.paco2 / 2, rel=0.15)


def test_tidal_volume_moves_plateau_pressure():
    p = pats()
    small = solve_steady(p, Vent(vt_ml_kg=4, rr=14))
    big = solve_steady(p, Vent(vt_ml_kg=10, rr=14))
    assert big.pplat > small.pplat + 5
    assert big.dp_stat > small.dp_stat


def test_peep_recruits_in_ards_and_the_cost_shows_up():
    p = ards()
    rows = [solve_steady(p, Vent(peep=peep, fio2=0.8, vt_ml_kg=6, rr=21))
            for peep in (5, 8, 10, 12, 15, 18, 22)]
    shunts = [r.shunt_pct for r in rows]
    assert shunts[0] > 15 and shunts[-1] < 8
    assert shunts[-1] < shunts[0]
    # recruitment is not a straight line: overdistension can push blood back
    # towards units that had just opened
    assert min(shunts[2:5]) <= shunts[1]
    stresses = [r.local_stress for r in rows]
    assert stresses[-1] > stresses[0]                              # overdistension
    cos = [r.co_l_min for r in rows]
    assert cos[-1] < cos[0]                                        # flow falls


def test_oxygenation_improves_with_fio2_and_shunt_blunts_it():
    healthy = pats()
    a = solve_steady(healthy, Vent(fio2=0.21, rr=12))
    b = solve_steady(healthy, Vent(fio2=0.5, rr=12))
    c = solve_steady(healthy, Vent(fio2=1.0, rr=12))
    assert a.pao2 < b.pao2 < c.pao2
    assert a.pao2 > 60                             # a normal lung still oxygenates
    sick = Patient(height_cm=170, weight_kg=70, severity=0.9, compliance=20,
                   heterogeneity=0.6)
    s_lo = solve_steady(sick, Vent(fio2=0.21, peep=5))
    s_hi = solve_steady(sick, Vent(fio2=1.0, peep=5))
    assert s_hi.pao2 > s_lo.pao2
    assert s_lo.pao2 < 60                          # the sick lung does not


def test_dynamic_hyperinflation_in_obstructive_disease():
    asthmatic = Patient(name="asthma", height_cm=175, weight_kg=80, severity=0.35,
                        compliance=55, resistance=22, exp_limit=2.5)
    tight = solve_steady(asthmatic, Vent(rr=24, ie=2.0, peep=5))
    slow = solve_steady(asthmatic, Vent(rr=24, ie=5.0, peep=5))
    assert tight.auto_peep > 5
    assert slow.auto_peep < tight.auto_peep
    assert tight.pplat > 30
    assert tight.co_l_min < slow.co_l_min          # trapping impedes venous return


def test_trapping_falls_as_exhalation_is_allowed_more_time():
    asthmatic = Patient(name="asthma", height_cm=175, weight_kg=80, severity=0.35,
                        compliance=55, resistance=22, exp_limit=2.5)
    values = [solve_steady(asthmatic, Vent(rr=20, ie=ie)).auto_peep
              for ie in (1.0, 1.5, 2.0, 3.0, 5.0)]
    assert all(a >= b - 1e-9 for a, b in zip(values, values[1:]))
    assert values[0] > values[-1]


def test_patent_foramen_ovale_bypasses_the_ventilator():
    p = pats(pfo=0.12)
    lo = solve_steady(p, Vent(peep=5, fio2=0.5))
    hi = solve_steady(p, Vent(peep=25, fio2=0.5))
    assert lo.shunt_pct > 10
    assert hi.shunt_pct >= lo.shunt_pct - 1e-9     # opening lung cannot fix a PFO
    assert hi.pao2 <= lo.pao2 + 1.0                # ... and barely helps


def test_capillary_loss_reads_as_dead_space_and_raises_co2():
    base = dict(name="pe", height_cm=170, weight_kg=70, severity=0.4, compliance=45,
                heterogeneity=0.6, perfusion_gradient=0.8)
    perfused = solve_steady(Patient(**base), Vent(peep=5, fio2=0.5, rr=16))
    blocked = solve_steady(Patient(**base, vascular_loss=0.6),
                           Vent(peep=5, fio2=0.5, rr=16))
    assert blocked.vd_vt > perfused.vd_vt + 0.10
    assert blocked.vd_vt > 0.40
    assert blocked.paco2 > perfused.paco2


def test_local_stress_exceeds_the_average_in_a_heterogeneous_lung():
    mild = ards(heterogeneity=0.0)
    rough = ards(heterogeneity=0.6)
    v = Vent(peep=10, vt_ml_kg=6, fio2=0.8, rr=21)
    rm, rr_ = solve_steady(mild, v), solve_steady(rough, v)
    # the worst unit feels far more than the average distending pressure, and
    # heterogeneity - not the average - is what pushes it up
    assert rm.local_stress > rm.dp_stat
    assert rr_.local_stress > rm.local_stress + 3.0
    assert rr_.pplat == pytest.approx(rm.pplat, rel=0.05)


def test_pressure_support_has_one_equilibrium_whatever_the_seed():
    p = Patient(name="awake", height_cm=172, weight_kg=72, severity=0.35,
                compliance=45, resistance=7, sedation=0.3, vco2=230.0)
    v = Vent(mode="psv", dp=12, rr=8, peep=8, fio2=0.4)
    a = solve_steady(p, v, paco2_state=20.0, sao2_state=0.99)
    b = solve_steady(p, v, paco2_state=80.0, sao2_state=0.90)
    assert a.paco2 == pytest.approx(b.paco2, abs=1e-3)
    assert a.drive == pytest.approx(b.drive, abs=1e-3)
    assert 10.0 < a.rr < 40.0 and a.effort_ratio > 0.9


def test_pressure_support_sedation_curve():
    def paco2_for(sed: float) -> float:
        p = Patient(name="awake", height_cm=172, weight_kg=72, severity=0.35,
                    compliance=45, resistance=7, sedation=sed, vco2=230.0)
        return solve_steady(p, Vent(mode="psv", dp=12, rr=8, peep=8, fio2=0.4)).paco2

    values = [paco2_for(s) for s in (0.0, 0.3, 0.6, 0.9, 1.0)]
    assert all(a <= b + 1e-6 for a, b in zip(values, values[1:]))
    assert values[0] < 45 and values[-1] > 70


def test_under_support_shows_as_effort():
    p = Patient(name="day4", height_cm=172, weight_kg=72, severity=0.30,
                compliance=45, resistance=7, sedation=0.15, muscle=0.8, vco2=240)
    weak = solve_steady(p, Vent(mode="psv", dp=4, rr=8, peep=5, fio2=0.35))
    enough = solve_steady(p, Vent(mode="psv", dp=12, rr=8, peep=5, fio2=0.35))
    assert weak.effort_ratio < 0.95 and weak.drive > 1.3
    assert enough.effort_ratio > 0.99 and enough.drive < weak.drive
    levels = {lvl for lvl, _ in result_flags(weak, p)}
    assert "alert" in levels or "warn" in levels


def test_clamping_keeps_nonsense_out():
    v = Vent(fio2=2.0, peep=-5, rr=400, vt_ml_kg=40, ie=0.1).clamped()
    assert 0.21 <= v.fio2 <= 1.0
    assert v.peep >= 0.0
    assert v.rr <= 60 and v.vt_ml_kg <= 15.0
    assert v.ie >= 0.35


def test_flags_fire_on_ventilator_harm():
    p = ards()
    r = solve_steady(p, Vent(vt_ml_kg=11, peep=10, fio2=0.8, rr=24))
    text = " ".join(t for _, t in r.flags)
    assert "driving pressure" in text or "Pplat" in text


def test_lung_units_span_apex_to_dependent():
    lung = Lung(n_units=12, severity=0.5)
    assert lung.z[0] < 0.1 and lung.z[-1] > 0.9
    assert sum(lung.perf) == pytest.approx(1.0, abs=0.35)
    assert any(lung.aff[i] > 0 for i in range(len(lung.z)))
