"""Generate reference cases (JSON) that any reimplementation must reproduce.

Used to check the JavaScript port in web/ against the Python model.  Both
steady-state operating points and time-stepped runs are emitted, because the
time course is the whole point of the tool.

Run:  python tools/gen_cases.py [outfile]
"""

from __future__ import annotations

import dataclasses
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ventsim import Patient, Simulation, Vent, solve_steady  # noqa: E402
from ventsim.scenarios import SCENARIOS  # noqa: E402

FIELDS = ("ph", "paco2", "pao2", "sao2", "hco3", "be", "pvco2", "pvo2", "vt_ml",
          "vt_ml_kg", "rr", "ti_s", "ve_l_min", "va_l_min", "vd_vt", "peco2",
          "peep", "auto_peep", "peep_total", "pplat", "pip", "dp_stat", "pmean",
          "peak_flow", "c_dyn", "local_stress", "mp_j_min", "fio2", "pf", "oi",
          "shunt_pct", "vq_low_pct", "vq_high_pct", "recruit_pct", "aa_gradient",
          "cao2", "cvao2_diff", "do2", "co_l_min", "drive", "effort_ratio")

SKIP = {"lung"}


def patient_dict(p: Patient) -> dict:
    out = {}
    for f in dataclasses.fields(p):
        if f.name in SKIP:
            continue
        v = getattr(p, f.name)
        if v is None and f.name not in ("vco2", "frc_l", "ve_rest"):
            continue
        out[f.name] = v
    # the "None" derived quantities are pinned to their resolved values so the
    # port does not have to guess how defaults are derived
    out["vco2"] = p.vco2
    out["frc_l"] = p.frc_l
    out["ve_rest"] = p.ve_rest
    return out


def vent_dict(v: Vent) -> dict:
    return {f.name: getattr(v, f.name) for f in dataclasses.fields(v)}


def expect(r) -> dict:
    return {k: getattr(r, k) for k in FIELDS}


cases: list[dict] = []


def steady(cid: str, p: Patient, v: Vent, hco3: float | None = None) -> None:
    kw = {} if hco3 is None else {"hco3": hco3}
    cases.append({"id": cid, "kind": "steady", "patient": patient_dict(p),
                  "vent": vent_dict(v), "hco3": hco3,
                  "expect": expect(solve_steady(p, v, **kw))})


def dynamic(cid: str, p: Patient, v: Vent, steps: list[tuple[float, dict]],
            dt: float = 0.05) -> None:
    """steps: (minutes_from_start, vent_changes) - sampled after each run.

    Contract for a reimplementation: apply `changes`, then advance exactly
    `round(span / dt)` steps of size `dt`, where `span` is the gap to the
    previous entry, and sample once at the end of each entry.
    """
    sim = Simulation(p, v)
    samples = []
    t_prev = 0.0
    for minutes, changes in steps:
        if changes:
            sim.set_vent(**changes)
        span = minutes - t_prev
        if span > 0:
            sim.run(span, sample_every=max(span, 1e-3), dt=dt)
        t_prev = minutes
        samples.append({"t_min": minutes, "changes": changes,
                        "expect": expect(sim.current())})
    cases.append({"id": cid, "kind": "dynamic", "patient": patient_dict(p),
                  "vent": vent_dict(v), "dt": dt, "start": expect(sim.current()),
                  "steps": samples})


# ---- every scenario at its own settings -----------------------------------
for s in SCENARIOS:
    sim = s.build()
    cases.append({"id": s.key, "kind": "steady", "patient": patient_dict(s.patient),
                  "vent": vent_dict(s.vent), "hco3": sim.hco3,
                  "expect": expect(sim.current())})

# ---- targeted sweeps ------------------------------------------------------
h = Patient(name="healthy", height_cm=172, weight_kg=75, compliance=60, resistance=5)
for fio2 in (0.21, 0.35, 1.0):
    steady(f"healthy_fio2_{fio2}", h, Vent(mode="vc", vt_ml_kg=8, rr=12, peep=5,
                                          fio2=fio2))
for vtk in (4, 6, 10, 14):
    steady(f"healthy_vt_{vtk}", h, Vent(mode="vc", vt_ml_kg=vtk, rr=12, peep=5, fio2=0.3))
for rr in (6, 30):
    steady(f"healthy_rr_{rr}", h, Vent(mode="vc", vt_ml_kg=8, rr=rr, peep=5, fio2=0.3))
steady("healthy_pc", h, Vent(mode="pc", dp=12, rr=14, ie=2.0, peep=5, fio2=0.4))
steady("healthy_psv", Patient(name="awake healthy", height_cm=172, weight_kg=75,
                              compliance=60, resistance=5, sedation=0.0),
       Vent(mode="psv", dp=10, rr=8, ie=2.0, peep=5, fio2=0.3))

ards = Patient(name="ARDS", height_cm=170, weight_kg=70, severity=0.55, compliance=28,
               resistance=8, heterogeneity=0.45, hgb=11.0, co_l_min=6.0,
               rv_sensitivity=1.2)
for peep in (5, 10, 15, 20, 25):
    steady(f"ards_peep_{peep}", ards,
           Vent(mode="vc", vt_ml_kg=6, rr=21, peep=peep, fio2=0.8))
for vtk in (4, 8):
    steady(f"ards_vt_{vtk}", ards,
           Vent(mode="vc", vt_ml_kg=vtk, rr=24, peep=12, fio2=0.7))
steady("ards_hb7", Patient(name="ARDS anaemic", height_cm=170, weight_kg=70,
                          severity=0.55, compliance=28, resistance=8, heterogeneity=0.45,
                          hgb=7.0, co_l_min=6.0, rv_sensitivity=1.2),
       Vent(mode="vc", vt_ml_kg=6, rr=21, peep=12, fio2=0.7))

copd = Patient(name="COPD", height_cm=170, weight_kg=70, severity=0.15, compliance=45,
               resistance=16, exp_limit=3.0, heterogeneity=0.45, vascular_loss=0.35,
               co_l_min=6.5, hgb=14.0, vco2=190.0)
for rr, ie in ((10, 4.0), (20, 1.5)):
    steady(f"copd_rr{rr}_ie{ie}", copd,
           Vent(mode="vc", vt_ml_kg=8, rr=rr, ie=ie, peep=5, fio2=0.3))

asthma = Patient(name="asthma", height_cm=178, weight_kg=75, severity=0.10,
                 compliance=55, resistance=30, exp_limit=6.0, heterogeneity=0.50,
                 co_l_min=7.0, vco2=210.0)
for rr, ie in ((12, 2.0), (8, 4.0)):
    steady(f"asthma_rr{rr}_ie{ie}", asthma,
           Vent(mode="vc", vt_ml_kg=8, rr=rr, ie=ie, peep=5, fio2=0.35))

pe = Patient(name="PE", height_cm=175, weight_kg=80, severity=0.05, compliance=55,
             resistance=6, heterogeneity=0.60, perfusion_gradient=0.8,
             vascular_loss=0.60, co_l_min=4.0, rv_sensitivity=2.0, vco2=230, pfo=0.06)
for peep in (5, 15):
    steady(f"pe_peep_{peep}", pe, Vent(mode="vc", vt_ml_kg=8, rr=16, peep=peep, fio2=0.4))

psv_pt = Patient(name="awake", height_cm=172, weight_kg=72, severity=0.30, compliance=45,
                 resistance=7, sedation=0.15, muscle=0.8, vco2=240, hgb=11.0, co_l_min=5.5)
for sed, dp in ((0.0, 8.0), (0.0, 20.0), (0.6, 12.0), (0.95, 12.0)):
    q = Patient(**{**patient_dict(psv_pt), "sedation": sed})
    steady(f"psv_sed{sed}_dp{dp}", q, Vent(mode="psv", dp=dp, rr=8, ie=2.0, peep=5, fio2=0.35))

steady("acidosis", Patient(name="DKA", height_cm=175, weight_kg=70, compliance=60,
                           resistance=5, hco3=12.0, acid_mmol_h=0.0),
       Vent(mode="vc", vt_ml_kg=8, rr=18, peep=5, fio2=0.3), hco3=12.0)
steady("renal_failure", Patient(name="CKD", height_cm=175, weight_kg=70, compliance=60,
                                resistance=5, renal=0.0, hco3=18.0),
       Vent(mode="vc", vt_ml_kg=8, rr=12, peep=5, fio2=0.3), hco3=18.0)

# ---- time courses ---------------------------------------------------------
hypo = Patient(name="hypercapnic", height_cm=172, weight_kg=75, severity=0.2,
               compliance=50, resistance=9, vco2=260.0)
dynamic("course_rr12to20", hypo, Vent(mode="vc", vt_ml_kg=8, rr=12, peep=6, fio2=0.4),
        [(0.5, {"rr": 20}), (2.0, {}), (10.0, {}), (40.0, {})])
dynamic("course_peep", ards, Vent(mode="vc", vt_ml_kg=6, rr=21, peep=8, fio2=0.8),
        [(0.0, {"peep": 18}), (1.0, {}), (5.0, {}), (30.0, {})], dt=0.1)
dynamic("course_fio2", pe, Vent(mode="vc", vt_ml_kg=8, rr=16, peep=5, fio2=0.3),
        [(0.0, {"fio2": 1.0}), (2.0, {}), (15.0, {})], dt=0.05)

chron = Patient(name="chronic", height_cm=168, weight_kg=65, severity=0.25,
                compliance=42, resistance=14, exp_limit=2.5, vco2=200)
chron_v = Vent(mode="vc", vt_ml_kg=6, rr=13, ie=3.0, peep=6, fio2=0.35)
sim = Simulation(chron, chron_v)
sim.adapt(48)
steady("chronic_48h_compensated", chron, chron_v, hco3=sim.hco3)
dynamic("course_post_hypercapnic", chron, chron_v,
        [(0.0, {"rr": 32}), (30.0, {}), (720.0, {})], dt=2.0)
dynamic("course_bicarb_load", Patient(**{**patient_dict(chron),
                                         "renal": 0.0, "acid_mmol_h": 8.0}),
        Vent(mode="vc", vt_ml_kg=8, rr=16, peep=5, fio2=0.3),
        [(30.0, {}), (180.0, {}), (720.0, {})], dt=2.0)

# ---- spontaneous breathing, stepped far coarser than the drive loop --------
# These catch a port that solves the drive/gas feedback in steps too large for
# it: the loop then limit-cycles instead of settling, and only these cases see
# it happen.
awake = Patient(name="awake", height_cm=172, weight_kg=75, sedation=0.0)
dynamic("course_psv_drive", awake,
        Vent(mode="psv", dp=12, rr=8, ie=2.0, peep=5, fio2=0.4),
        [(5.0, {}), (15.0, {})], dt=5.0)
dynamic("course_psv_step_down", psv_pt,
        Vent(mode="psv", dp=14, rr=8, ie=2.0, peep=5, fio2=0.35),
        [(0.0, {"dp": 6.0}), (10.0, {})], dt=1.0)

out = Path(sys.argv[1] if len(sys.argv) > 1 else "web/cases.json")
out.parent.mkdir(parents=True, exist_ok=True)
# write-then-rename so a concurrent reader never sees a half-written file
tmp = out.with_suffix(out.suffix + ".tmp")
tmp.write_text(json.dumps({"version": 1, "n_cases": len(cases), "cases": cases}) + "\n")
os.replace(tmp, out)
print(f"{len(cases)} cases -> {out}")
