#!/usr/bin/env python3
"""Compare a scenario's JS Simulation state with Python's, field by field.

    python3 tools/probe_js.py chronic_hypercapnic
"""
from __future__ import annotations

import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ventsim.scenarios import get  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
URL = (ROOT / "web" / "ventsim.html").as_uri()
KEY = sys.argv[1] if len(sys.argv) > 1 else "chronic_hypercapnic"

FIELDS = ["hco3", "drive"]
RESULT_FIELDS = ["ph", "paco2", "pao2", "sao2", "hco3", "pplat", "dp_stat", "pip",
                "pmean", "auto_peep", "recruit_pct", "shunt_pct", "co_l_min",
                "vt_ml", "ve_l_min", "drive", "effort_ratio", "local_stress",
                "stress_pct", "mp_j_min", "vd_vt", "pvco2", "pvo2", "cao2", "cvao2_diff"]

PROBE = """(key) => {
  loadScenario(key);
  const out = {hco3: sim.hco3, drive: sim.drive, t: sim.t_min,
               fastCO2: sim._fastCO2 !== undefined ? sim._fastCO2 : sim.fastCO2,
               paco2: sim.paco2, pao2: sim.pao2};
  const tgt = sim.target(), cur = sim.current();
  const res = {};
  const fields = %s;
  for (const f of fields) { res['target.' + f] = tgt[f]; res['current.' + f] = cur[f]; }
  const p = sim.patient;
  for (const f of ['hco3','vco2','severity','heterogeneity','perfusion_gradient',
                   'vascular_loss','n_units','pfo','base_shunt','compliance','resistance',
                   'sedation','set_point','muscle','exp_limit','affected_fraction'])
    res['patient.' + f] = p[f];
  const v = sim.vent;
  for (const f of ['mode','fio2','peep','vt_ml','rr','ie','dp']) res['vent.' + f] = v[f];
  return {sim: out, res};
}""" % (str(RESULT_FIELDS),)


def main() -> int:
    scen = get(KEY)
    sim = scen.build()
    tgt, cur = sim.target(), sim.current()
    with sync_playwright() as pw:
        b = pw.chromium.launch(channel="chrome", args=["--no-sandbox", "--disable-gpu"])
        pg = b.new_page()
        pg.goto(URL)
        pg.wait_for_timeout(300)
        got = pg.evaluate(PROBE, KEY)
        b.close()

    print(f"scenario {KEY}: js.sim / py.sim")
    print(f"  hco3   {got['sim']['hco3']!r:>22} / {sim.hco3!r}")
    print(f"  drive  {got['sim']['drive']!r:>22} / {sim.drive!r}")
    print(f"  paco2  {got['sim']['paco2']!r:>22} / {sim.paco2!r}")
    print(f"  pao2   {got['sim']['pao2']!r:>22} / {sim.pao2!r}")
    print(f"  t      {got['sim']['t']!r:>22} / {sim.t_min!r}")
    worst = 0.0
    for f in RESULT_FIELDS:
        jv = got["res"].get("target." + f)
        pv = getattr(tgt, f, None)
        cv = getattr(cur, f, None)
        d = abs(jv - pv) if isinstance(jv, (int, float)) else float("nan")
        dc = abs(jv - cv) if isinstance(jv, (int, float)) else float("nan")
        worst = max(worst, d)
        flag = "" if d < 1e-6 else ("   <-- target differs" if d > 1e-3 else "")
        print(f"  {f:<14} {jv!r:>22} / target {pv!r}   (|Δt|={d:.6g}) cur {cv!r} |Δc|={dc:.6g}{flag}")
    print("\npatient/vent resolution:")
    for k, v in sorted(got["res"].items()):
        if k.startswith("patient.") or k.startswith("vent."):
            print(f"  {k} = {v!r}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
