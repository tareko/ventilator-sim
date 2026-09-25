"""Time course: how long the changes take, and what the slow processes do."""

import pytest

from ventsim import Patient, Simulation, Vent


def hypercapnic(**kw) -> Simulation:
    p = Patient(name="hypercapnic", height_cm=170, weight_kg=70, severity=0.3,
                compliance=40, resistance=9, vco2=230.0, **kw)
    sim = Simulation(p, Vent(rr=8, vt_ml_kg=7, peep=5, fio2=0.4))
    for _ in range(200):
        sim.step(0.5)
    return sim


def test_the_lung_blood_store_fills_within_minutes():
    sim = hypercapnic()
    start = sim.current().paco2
    assert start > 55                                    # CO2 has accumulated
    sim.set_vent(rr=20)
    tgt = sim.target().paco2
    assert tgt < 40
    seen = [start]
    for _ in range(8):                                    # two minutes
        sim.step(0.25)
        seen.append(sim.current().paco2)
    fast_phase = seen[0] - seen[-1]
    for _ in range(13):                                   # another quarter hour
        sim.step(1.0)
        seen.append(sim.current().paco2)
    assert all(a >= b for a, b in zip(seen, seen[1:]))     # monotone fall
    assert fast_phase > 0.5 * (seen[0] - seen[-1])         # most of it at once
    # the last couple of mmHg belong to the tissue store and take longer
    assert seen[-1] == pytest.approx(tgt, abs=3.0)


def test_the_tissue_store_makes_it_slow_at_the_end():
    sim = hypercapnic()
    sim.set_vent(rr=20)
    after_2 = sim.run(2, sample_every=2, dt=0.1)[-1].paco2
    after_30 = sim.run(28, sample_every=2, dt=0.1)[-1].paco2
    tgt = sim.target().paco2
    assert after_2 > tgt + 1.0                            # not there yet
    assert after_30 == pytest.approx(tgt, abs=1.0)         # ... getting there
    assert after_2 > after_30


def test_hco3_property_mixes_set_point_and_acute_buffering():
    sim = Simulation(Patient(hco3=24.0), Vent(rr=6, vt_ml_kg=7))
    assert sim.hco3 == pytest.approx(24.0, abs=0.01)
    sim.run(10, sample_every=10, dt=0.5)
    with_co2 = sim.hco3
    # CO2 is high, so cellular buffering has added bicarbonate on top of the
    # metabolic set point, which has not moved
    assert sim.paco2 > 55
    assert with_co2 > 24.0
    assert sim._hco3_base == pytest.approx(24.0, abs=0.4)
    sim.hco3 = 30.0                                        # the set point moves
    assert sim.hco3 == pytest.approx(30.0, abs=1e-9)


def test_renal_compensation_needs_perfusion_and_time():
    sim = hypercapnic()
    before = sim.hco3
    samples = sim.adapt(48)
    assert len(samples) == 4
    assert sim.hco3 > before + 3.0
    assert sim.current().ph > 7.32                          # still acidotic, less so

    failing = hypercapnic(renal=0.0)
    failing.adapt(48)
    assert failing.hco3 < before + 1.5                      # no kidneys, no compensation


def test_post_hypercapnic_alkalosis():
    sim = hypercapnic()
    sim.adapt(48)
    compensated = sim.hco3
    sim.set_vent(rr=26)
    sim.run(20, sample_every=10, dt=0.2)
    r = sim.current()
    assert r.paco2 < 35
    assert r.hco3 > 30 and r.hco3 < compensated             # buffering went, kidneys did not
    assert r.ph > 7.60                                      # dangerous after a long stay


def test_acid_and_bicarbonate_loads_move_the_base_line():
    p = Patient(name="dka", height_cm=170, weight_kg=70, acid_mmol_h=20.0,
                sedation=1.0)
    sim = Simulation(p, Vent(rr=12, vt_ml_kg=8))
    sim.run(12 * 60, sample_every=60, dt=2.0)
    acidotic = sim.current()
    assert acidotic.hco3 < 21 and acidotic.ph < 7.35
    assert acidotic.be < -4

    loaded = Simulation(Patient(bicarb_mmol_h=10.0, sedation=1.0), Vent(rr=14, vt_ml_kg=8))
    loaded.run(12 * 60, sample_every=60, dt=2.0)
    assert loaded.current().hco3 > 25
    assert loaded.current().be > 2


def test_drive_relaxes_toward_the_response_line():
    p = Patient(name="awake", height_cm=172, weight_kg=75, sedation=0.2)
    sim = Simulation(p, Vent(mode="psv", dp=14, rr=8, peep=5))
    sim.run(30, sample_every=10, dt=0.2)
    assert sim.drive == pytest.approx(1.0, abs=0.25)        # near its resting minute ventilation
    # and a deep sedation removes it
    sed = Simulation(Patient(sedation=1.0), Vent(mode="psv", dp=14, rr=8))
    sed.run(30, sample_every=10, dt=0.2)
    assert sed.drive == pytest.approx(0.0, abs=1e-6)


def test_run_samples_and_clock():
    sim = Simulation(Patient(), Vent())
    assert sim.run(0) == []
    rows = sim.run(10, sample_every=2, dt=0.1)
    assert len(rows) == 5
    assert sim.t_min == pytest.approx(10.0, abs=1e-6)
    assert rows[-1].t_min == pytest.approx(sim.t_min, abs=1e-9)


def test_step_is_stable_for_a_coarse_interval():
    sim = Simulation(Patient(sedation=0.0), Vent(mode="psv", dp=10, rr=8))
    for _ in range(40):
        sim.step(5.0)
        r = sim.current()
        assert 4 < r.paco2 < 160
        assert 7.0 < r.ph < 7.60
        assert r.sao2 >= 0.0 and r.sao2 <= 1.0
