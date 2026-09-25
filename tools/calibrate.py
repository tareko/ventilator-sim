"""Calibration / sanity battery: does the model behave like a bedside teaching case?

Run with:  python tools/calibrate.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ventsim import Patient, Simulation, Vent, solve_steady  # noqa: E402


def line(tag: str, r) -> str:
    return (f"{tag:<26} pH {r.ph:5.2f}  PaCO2 {r.paco2:5.1f}  PaO2 {r.pao2:6.1f}  "
            f"SaO2 {100*r.sao2:5.1f}%  HCO3 {r.hco3:4.1f}   "
            f"Pplat {r.pplat:5.1f} dP {r.dp_stat:4.1f} autoPEEP {r.auto_peep:4.1f} "
            f"Pmean {r.pmean:5.1f}  Vd/Vt {r.vd_vt:4.2f} shunt {r.shunt_pct:4.1f}% "
            f"CO {r.co_l_min:3.1f}  A-a {r.aa_gradient:5.1f}")


print("=" * 118)
print("1. healthy adult, sedated, volume control, room air  (want pH ~7.40, PaCO2 ~40, "
      "PaO2 ~95, Pplat ~15)")
print("=" * 118)
p0 = Patient(name="healthy", height_cm=172, weight_kg=75, compliance=60, resistance=5)
r = solve_steady(p0, Vent(mode="vc", vt_ml_kg=8, rr=12, peep=5, fio2=0.21))
print(line("VC 8ml/kg RR12", r))
r = solve_steady(p0, Vent(mode="vc", vt_ml_kg=8, rr=12, peep=5, fio2=1.0))
print(line("same, FiO2 1.0", r))
for rr in (8, 10, 12, 16, 20):
    r = solve_steady(p0, Vent(mode="vc", vt_ml_kg=8, rr=rr, peep=5, fio2=0.21))
    print(line(f"RR {rr}", r))
for vtk in (4, 6, 8, 10, 12):
    r = solve_steady(p0, Vent(mode="vc", vt_ml_kg=vtk, rr=12, peep=5, fio2=0.21))
    print(line(f"Vt {vtk} ml/kg", r))

print()
print("=" * 118)
print("2. moderate ARDS: PEEP sweep at FiO2 0.9  (want a PaO2 peak, rising dead space "
      "and falling CO high up)")
print("=" * 118)
ards = Patient(name="ARDS mod", severity=0.55, compliance=28, resistance=8,
               co_l_min=6.0, rv_sensitivity=1.2, hgb=11.0, weight_kg=70, height_cm=170)
for peep in (5, 8, 10, 12, 15, 18, 22, 26, 30):
    r = solve_steady(ards, Vent(mode="vc", vt_ml_kg=6, rr=28, peep=peep, fio2=0.9))
    print(f"PEEP {peep:2d}  " + line("", r))

print()
print("=" * 118)
print("3. Vt sweep in ARDS at PEEP 15 (the 'baby lung': local stress climbs fast)")
print("=" * 118)
for vtk in (4, 5, 6, 7, 8, 10):
    r = solve_steady(ards, Vent(mode="vc", vt_ml_kg=vtk, rr=28, peep=15, fio2=0.9))
    print(f"Vt {vtk:2d}   " + line("", r))

print()
print("=" * 118)
print("4. COPD / emphysema: rate and I:E (dynamic hyperinflation)")
print("=" * 118)
copd = Patient(name="COPD", severity=0.15, compliance=45, resistance=16, exp_limit=3.0,
               hgb=14.0, co_l_min=6.5, vco2=180.0, hco3=24.0, sedation=1.0,
               heterogeneity=0.45, vascular_loss=0.35)
for rr, ie in ((12, 2.0), (12, 3.0), (10, 3.5), (14, 2.0), (8, 3.0), (20, 2.0)):
    r = solve_steady(copd, Vent(mode="vc", vt_ml_kg=8, rr=rr, ie=ie, peep=5, fio2=0.30))
    print(f"RR {rr:2d} I:E 1:{ie:g}  " + line("", r))

print()
print("=" * 118)
print("5. severe asthma: the same settings destroy you")
print("=" * 118)
asthma = Patient(name="asthma", severity=0.10, compliance=55, resistance=30, exp_limit=6.0,
                 hgb=15, co_l_min=7.0, vco2=210.0, heterogeneity=0.5)
for rr, ie in ((12, 2.0), (10, 4.0), (8, 4.0)):
    r = solve_steady(asthma, Vent(mode="vc", vt_ml_kg=8, rr=rr, ie=ie, peep=5, fio2=0.35))
    print(f"RR {rr:2d} I:E 1:{ie:g}  " + line("", r))

print()
print("=" * 118)
print("6. PE: dead-space dominant, PEEP does nothing good")
print("=" * 118)
pe = Patient(name="PE", severity=0.05, compliance=55, resistance=6, heterogeneity=0.6,
             perfusion_gradient=0.8, vascular_loss=0.6, co_l_min=4.0, rv_sensitivity=2.0,
             vco2=230)
for peep in (5, 10, 15):
    r = solve_steady(pe, Vent(mode="vc", vt_ml_kg=8, rr=16, peep=peep, fio2=0.6))
    print(f"PEEP {peep:2d}  " + line("", r))

print()
print("=" * 118)
print("7. spontaneous breathing (PSV): sedation vs CO2")
print("=" * 118)
sb = Patient(name="awake", severity=0.35, compliance=45, resistance=7, sedation=0.0,
             vco2=230.0)
for sed in (0.0, 0.3, 0.6, 0.9, 1.0):
    q = Patient(name="awake", severity=0.35, compliance=45, resistance=7, sedation=sed,
                vco2=230.0)
    r = solve_steady(q, Vent(mode="psv", dp=12, rr=8, peep=8, fio2=0.4))
    print(f"sedation {sed:.1f} " + line("", r))

print()
print("=" * 118)
print("8. time course: raise RR 12->20 in a hypercapnic patient (watch PaCO2 fall)")
print("=" * 118)
hypo = Patient(name="hypercapnic", severity=0.2, compliance=50, resistance=9,
               vco2=260.0, sedation=1.0)
s = Simulation(hypo, Vent(mode="vc", vt_ml_kg=8, rr=12, peep=6, fio2=0.4))
print(f"start      " + line("", s.current()))
s.set_vent(rr=20)
last = s.t_min
for target in (0.5, 1, 2, 5, 10, 20, 40):
    s.run(target - last, sample_every=max(target - last, 0.05), dt=0.05)
    last = target
    print(f"t={target:4.1f} min " + line("", s.current()))

print()
print("=" * 118)
print("9. chronic compensation: 48 h at PaCO2 ~60, then 'normalise' the CO2 to 40")
print("=" * 118)
chron = Patient(name="chronic hypercapnic", severity=0.25, compliance=42, resistance=14,
                exp_limit=2.5, vco2=200.0, sedation=1.0)
s2 = Simulation(chron, Vent(mode="vc", vt_ml_kg=6, rr=13, ie=3.0, peep=6, fio2=0.35))
print(f"baseline   " + line("", s2.current()))
s2.adapt(48)
print(f"after 48 h " + line("", s2.current()) + f"   t={s2.t_min/60:.1f} h")
s2.set_vent(rr=32)
s2.run(20, sample_every=20, dt=0.2)
print(f"acute hypo " + line("", s2.current()))
s2.run(48 * 60, sample_every=24 * 60, dt=5)
print(f"+48 h later" + line("", s2.current()))
