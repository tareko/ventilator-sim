#!/usr/bin/env python3
"""Drive web/ventsim.html with Playwright and check it against the Python model.

Loads a scenario, applies a few setting changes through the real control
listeners, runs the clock, and prints what the DOM shows next to what Python
computes for the same sequence.

    python3 tools/check_ui.py [scenario_key]
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ventsim.scenarios import get  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
URL = (ROOT / "web" / "ventsim.html").as_uri()
KEY = sys.argv[1] if len(sys.argv) > 1 else "ards_mod"

JS_PANEL = """(sel) => {
  const e = document.querySelector(sel);
  return e ? e.innerText : '';
}"""

JS_SET = """([group, field, value]) => {
  const ctl = controls[group][field];
  if (!ctl) return 'no control ' + group + '.' + field;
  ctl.number.value = String(value);
  ctl.number.dispatchEvent(new Event('input'));
  ctl.number.dispatchEvent(new Event('change'));
  return 'ok';
}"""

LABELS = {"ph": "ph", "paco2": "paco2", "pao2": "pao2", "sao2": "sao2",
          "hco3": "hco3", "plateau": "pplat", "driving pressure": "dp_stat",
          "peak pressure": "pip", "auto-peep": "auto_peep",
          "shunt %": "shunt_pct", "recruited %": "recruit_pct",
          "cardiac output l/min": "co_l_min", "mean airway": "pmean"}

# fields compared between the DOM and the Python model
want = ("ph", "paco2", "pao2", "sao2", "hco3", "pplat", "dp_stat")


def parse_panel(text: str) -> dict[str, float]:
    """Pull label/number pairs out of a rendered panel."""
    out: dict[str, float] = {}
    tokens = re.split(r"\s*\n\s*", text.strip())
    for i, tok in enumerate(tokens):
        key = tok.strip().lower().rstrip("%").strip()
        if key in LABELS and i + 1 < len(tokens):
            m = re.match(r"^(-?\d+(?:\.\d+)?)", tokens[i + 1].strip())
            if m:
                val = float(m.group(1))
                if key == "sao2":
                    val /= 100.0
                out[LABELS[key]] = val
    return out


def main() -> int:
    scen = get(KEY)
    sim = scen.build()
    worst = 0.0
    with sync_playwright() as pw:
        b = pw.chromium.launch(channel="chrome", args=["--no-sandbox", "--disable-gpu"])
        pg = b.new_page(viewport={"width": 1500, "height": 1150})
        problems: list[str] = []
        pg.on("pageerror", lambda e: problems.append(str(e)))
        pg.goto(URL)
        pg.wait_for_timeout(400)
        print(f"page: {pg.title()}")
        pg.select_option("#scenario", KEY)
        pg.wait_for_timeout(300)

        def show(tag: str, ref_result) -> dict[str, float]:
            abg = pg.evaluate(JS_PANEL, "#abg")
            der = pg.evaluate(JS_PANEL, "#derived")
            got = parse_panel(abg)
            got.update(parse_panel(der))
            ref = {k: getattr(ref_result, k) for k in want}
            deltas = {k: abs(got[k] - ref[k]) for k in got if k in want}
            bad = max(deltas.values()) if deltas else 0.0
            nonlocal worst
            worst = max(worst, bad)
            print(f"\n{tag}")
            print("  ui    :", " ".join(f"{k}={got[k]:.4g}" for k in sorted(got)))
            print("  python:", " ".join(f"{k}={ref[k]:.4f}" for k in sorted(ref)))
            print(f"  worst |delta| = {bad:.4f}")
            return got

        print("siteline:", pg.evaluate(JS_PANEL, "#siteline"))
        show(f"scenario {KEY} baseline (operating point)", sim.target())

        for group, field, value, label in [
            ("vent", "peep", 18, "PEEP 18"),
            ("vent", "fio2", 0.5, "FiO2 0.5"),
        ]:
            if pg.evaluate(JS_SET, [group, field, value]) != "ok":
                problems.append(f"control {group}.{field} missing")
                continue
            pg.wait_for_timeout(150)
            sim.set_vent(**{field: value})
            show(f"after {label} (operating point)", sim.target())

        # switch the panel to the simulated "now" state and run the clock
        try:
            pg.click("#viewseg >> text=now")
        except Exception:
            problems.append("could not switch panel to 'now' view")
        pg.click("#step1")
        pg.wait_for_timeout(200)
        sim.step(0.1)
        show("after one 0.1 min step (now)", sim.current())

        pg.fill("#mins", "10")
        pg.click("#run")
        pg.wait_for_timeout(400)
        sim.run(10, sample_every=10, dt=0.1)
        print("\nclock:", pg.evaluate(JS_PANEL, "#clock"))
        show("after running 10 min (now)", sim.current())
        print("trends:", pg.evaluate(JS_PANEL, "#trends")[:400])
        shots = ROOT / ".artifacts"
        shots.mkdir(exist_ok=True)
        pg.screenshot(path=str(shots / f"ui_{KEY}.png"), full_page=True)
        b.close()
        print("\nJS page errors:", problems or "none")
        print(f"\nWORST UI/Python delta across the session: {worst:.4f}")
        return 0 if worst < 0.05 and not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())
