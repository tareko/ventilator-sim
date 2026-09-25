"""Command-line bedside shell: change settings, read the blood gas, run the clock.

    python -m ventsim                                    # interactive
    python -m ventsim --scenario ards_mod --optimize     # ranked settings
    python -m ventsim --scenario copd --set rr=20 ie=3 --watch 30m

Type `help` inside the shell for the command list.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path

from . import __version__
from .dynamics import Simulation
from .model import Patient, Result, Vent, ibw_kg
from .optimizer import Weights, cost_of, explain, optimize
from .physiology import so2
from .scenarios import BY_KEY, SCENARIOS, Scenario

# --------------------------------------------------------------- parsing

_DURATION = re.compile(r"^(\d+(?:\.\d*)?)\s*([smh]?)$")
MODES = ("vc", "pc", "psv")


def parse_duration(text: str) -> float:
    """'90s', '10m', '2h', '1.5' -> minutes."""
    m = _DURATION.match(str(text).strip())
    if not m:
        raise ValueError(f"not a duration: {text!r} (use 45s, 10m, 2h)")
    val, unit = float(m.group(1)), (m.group(2) or "m")
    return val * {"s": 1.0 / 60.0, "m": 1.0, "h": 60.0}[unit]


# user-facing key -> attribute on Vent / Patient
VENT_KEYS = {
    "mode": "mode",
    "fio2": "fio2",
    "peep": "peep",
    "vt": "vt_ml_kg",
    "vt_ml_kg": "vt_ml_kg",
    "rr": "rr",
    "rate": "rr",
    "ie": "ie",
    "ie_ratio": "ie",
    "dp": "dp",
    "ps": "dp",
    "pressure_support": "dp",
    "ti": "ti_s",
    "ti_s": "ti_s",
}
PATIENT_KEYS = {
    "name": "name",
    "sex": "sex",
    "height": "height_cm",
    "weight": "weight_kg",
    "compliance": "compliance",
    "resistance": "resistance",
    "exp_limit": "exp_limit",
    "severity": "severity",
    "affected": "affected_fraction",
    "affected_fraction": "affected_fraction",
    "heterogeneity": "heterogeneity",
    "scatter": "heterogeneity",
    "perf_gradient": "perfusion_gradient",
    "perfusion_gradient": "perfusion_gradient",
    "vessels": "vascular_loss",
    "vascular_loss": "vascular_loss",
    "base_shunt": "base_shunt",
    "pfo": "pfo",
    "hpv": "hpv",
    "n_units": "n_units",
    "hgb": "hgb",
    "hb": "hgb",
    "co": "co_l_min",
    "cardiac_output": "co_l_min",
    "rv_sensitivity": "rv_sensitivity",
    "vco2": "vco2",
    "rq": "rq",
    "temp": "temp_c",
    "hco3": "hco3",
    "renal": "renal",
    "acid": "acid_mmol_h",
    "acid_mmol_h": "acid_mmol_h",
    "bicarb_load": "bicarb_mmol_h",
    "bicarb_mmol_h": "bicarb_mmol_h",
    "sedation": "sedation",
    "sed": "sedation",
    "set_point": "set_point",
    "muscle": "muscle",
    "ve_rest": "ve_rest",
    "rr_rest": "rr_rest",
    "notes": "notes",
}
INT_PATIENT = {"n_units"}
TEXT_PATIENT = {"name", "notes", "sex"}


@dataclass
class Assignment:
    key: str          # "vent.peep" or "patient.hgb"
    value: object


def parse_assignments(tokens: list[str], *, patient_only: bool = False) -> list[Assignment]:
    """Turn ['peep=12', 'sedation=0.5'] into typed assignments."""
    out: list[Assignment] = []
    for tok in tokens:
        if "=" not in tok:
            raise ValueError(f"expected key=value, got {tok!r}")
        key, _, raw = tok.partition("=")
        key = key.strip().lower()
        raw = raw.strip()
        attr_v = None if patient_only else VENT_KEYS.get(key)
        attr_p = PATIENT_KEYS.get(key)
        if attr_v is None and attr_p is None:
            raise ValueError(
                f"unknown setting {key!r}.  Ventilator: mode fio2 peep vt rr ie dp ti."
                "  Patient: compliance resistance severity scatter hgb co sedation"
                " muscle renal hco3 acid bicarb_load vco2 rq temp exp_limit"
                " perf_gradient vessels pfo hpv n_units height weight sex name notes")
        namespace, attr = ("vent", attr_v) if attr_v is not None else ("patient", attr_p)
        if attr == "mode":
            value: object = raw.lower()
            if value not in MODES:
                raise ValueError("mode is vc, pc or psv")
        elif attr in TEXT_PATIENT:
            value = raw
        else:
            try:
                num = float(raw)
            except ValueError:
                raise ValueError(f"{key} needs a number, got {raw!r}") from None
            if attr in INT_PATIENT:
                value = int(round(num))
            elif attr == "ti_s" and num <= 0:
                value = None                # ti=0 means "let I:E decide"
            else:
                value = num
        out.append(Assignment(f"{namespace}.{attr}", value))
    return out


KNOB_ALIASES = {"peep": "peep", "fio2": "fio2", "vt": "vt_ml_kg",
                "vt_ml_kg": "vt_ml_kg", "rr": "rr", "rate": "rr",
                "ie": "ie", "dp": "dp", "ps": "dp"}


def grid_from_tokens(tokens: list[str]) -> dict:
    """'vt=4,6' -> {'vt_ml_kg': [4.0, 6.0]}, with the shell's aliases."""
    grids: dict[str, list[float]] = {}
    for tok in tokens:
        key, sep, vals = tok.partition("=")
        attr = KNOB_ALIASES.get(key.strip().lower())
        if attr is None:
            raise ValueError(f"cannot optimise {key!r};"
                             " knobs are peep fio2 vt rr ie dp")
        if not sep:
            raise ValueError(f"expected {key}=v1,v2,...")
        try:
            nums = [float(x) for x in vals.split(",") if x.strip()]
        except ValueError:
            raise ValueError(f"bad values in {tok!r}") from None
        if not nums:
            raise ValueError(f"no values for {key!r}")
        grids.setdefault(attr, []).extend(nums)
    if not grids:
        raise ValueError("empty search grid")
    return grids


# --------------------------------------------------------------- rendering

def fmt(value: float | None, width: int = 6, digits: int = 1) -> str:
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return " " * (width - 1) + "-"
    return f"{value:>{width}.{digits}f}"


LEVEL_MARK = {"alert": "!", "warn": "~", "note": "."}


def ibw_of(patient: Patient) -> float:
    return ibw_kg(patient.height_cm, patient.sex)


def settings_line(vent: Vent, patient: Patient) -> str:
    ibw = ibw_of(patient)
    if vent.mode == "vc":
        mid = (f"Vt {vent.vt_ml_kg:g} ml/kg = {vent.vt_ml_kg * ibw:.0f} ml"
               f"   RR {vent.rr:g}   I:E 1:{vent.ie:g}")
    elif vent.mode == "pc":
        mid = f"dP {vent.dp:g}   RR {vent.rr:g}   I:E 1:{vent.ie:g}"
    else:
        mid = f"PS {vent.dp:g}   RR {vent.rr:g}   I:E 1:{vent.ie:g}"
    return (f"{vent.label}   FiO2 {vent.fio2:.2f}   PEEP {vent.peep:g}   "
            f"{mid}")


def patient_line(patient: Patient) -> str:
    p = patient
    return (f"{p.name}   IBW {ibw_of(p):.1f} kg   wt {p.weight_kg:g} kg   "
            f"C {p.compliance:g}   R {p.resistance:g}   severity {p.severity:g}"
            f"   Hb {p.hgb:g}   CO {p.co_l_min:g} L/min   sedation {p.sedation:g}")


def gas_block(r: Result) -> str:
    rows = [
        (f"pH {r.ph:.2f}   PaCO2 {fmt(r.paco2, 5)}   PaO2 {fmt(r.pao2, 6)}"
         f"   SaO2 {100 * r.sao2:5.1f}%   HCO3 {r.hco3:.1f}   BE {r.be:+.1f}"),
        (f"PvCO2 {r.pvco2:.1f}   PvO2 {r.pvo2:.1f}   A-a {r.aa_gradient:.1f}"
         f"   P/F {r.pf:.0f}   OI {r.oi:.1f}   shunt {r.shunt_pct:.1f}%"
         f"   Vd/Vt {r.vd_vt:.2f}"),
        (f"Pplat {r.pplat:.1f}   Ppeak {r.pip:.1f}   dP {r.dp_stat:.1f}"
         f"   Pmean {r.pmean:.1f}   auto-PEEP {r.auto_peep:.1f}"
         f"   Ti {r.ti_s:.2f} s"),
        (f"Vt {r.vt_ml:.0f} ml ({r.vt_ml_kg:.1f}/kg)   RR {r.rr:g}"
         f"   VE {r.ve_l_min:.2f} L/min   VA {r.va_l_min:.2f}"
         f"   MP {r.mp_j_min:.1f} J/min   Cdyn {r.c_dyn:.0f}"),
        (f"recruited {r.recruit_pct:.0f}%   low V/Q {r.vq_low_pct:.1f}%"
         f"   high V/Q {r.vq_high_pct:.1f}%   local stress {r.local_stress:.1f}"
         f"   CO {r.co_l_min:.1f} L/min   drive {r.drive:.2f}x"),
    ]
    return "\n".join("  " + row for row in rows)


def flags_block(r: Result) -> str:
    if not r.flags:
        return "  flags     nothing flagged"
    lines = []
    for i, (level, text) in enumerate(r.flags):
        head = "  flags     " if i == 0 else "              "
        lines.append(f"{head}{LEVEL_MARK.get(level, '.')} {text}")
    return "\n".join(lines)


def effort_note(r: Result) -> str:
    if r.mode != "psv":
        return ""
    return (f"  effort    delivering {100 * r.effort_ratio:.0f}% of the tidal"
            " volume the patient is asking for")


def trend_table(rows: list[Result]) -> str:
    out = ["      t       pH   PaCO2    PaO2   SaO2     HCO3    Pplat"
           "      dP  aPEEP      MP     CO  drive"]
    for r in rows:
        out.append(
            f"{r.t_min:>8.1f} {r.ph:>7.3f} {r.paco2:>7.1f} {r.pao2:>8.1f}"
            f" {100 * r.sao2:>6.1f} {r.hco3:>8.2f} {r.pplat:>8.1f}"
            f" {r.dp_stat:>7.1f} {r.auto_peep:>6.1f} {r.mp_j_min:>7.1f}"
            f" {r.co_l_min:>6.1f} {r.drive:>6.2f}")
    return "\n".join(out)


def vq_block(r: Result) -> str:
    """Per-unit table plus a perfusion histogram over V/Q bins."""
    units = list(r.units or [])
    if not units:
        return "  no per-unit data"
    co = max(r.co_l_min, 1e-9)
    q_tot_l = max(sum(u.perf for u in units) * co, 1e-9)
    va_max = max((u.va for u in units), default=1.0)
    p_max = max((u.perf for u in units), default=1.0)
    lines = ["  unit      z   perf    vent       V/Q   PcCO2   PcO2"
             "  sat   O2cont  open  over   vent (V) and perf (P)",
             "              L/min mL/min             mmHg  mmHg"
             "     mL/mL"]
    for i, u in enumerate(units, start=1):
        nv = int(round(18 * u.va / va_max)) if va_max > 0 else 0
        nq = int(round(18 * u.perf / p_max)) if p_max > 0 else 0
        bar = "".join("V" if j < nv and j >= nq else
                      "B" if j < nv and j < nq else
                      "P" if j < nq else " " for j in range(18))
        sat = 100.0 * so2(u.po2, None, r.ph, r.paco2)
        lines.append(
            f"{i:>4d} {u.z:>6.2f} {u.perf * co:>7.2f} {u.va:>7.1f}"
            f" {u.vq:>9.3f}"
            f" {u.pc_co2:>7.1f} {u.po2:>6.1f} {sat:>5.0f}% {u.cc_o2:>7.3f}"
            f" {100 * u.recruited:>5.0f}% {100 * u.overdist:>5.0f}%  [{bar}]")
    lines.append("")
    lines.append(f"  perfusion by V/Q bin (total perfusion {q_tot_l:.2f} L/min,"
                 f" {len(units)} units)")
    bins = [(0.0, 0.05, "shunt     "), (0.05, 0.30, "V/Q <0.3  "),
            (0.30, 0.60, "V/Q 0.3-0.6"), (0.60, 1.00, "V/Q 0.6-1 "),
            (1.00, 1.50, "V/Q 1-1.5 "), (1.50, 3.00, "V/Q 1.5-3 "),
            (3.00, float("inf"), "V/Q >3    ")]
    for lo, hi, label in bins:
        in_bin = [u for u in units if lo <= u.vq < hi]
        q = sum(u.perf for u in in_bin)
        frac = q / max(sum(u.perf for u in units), 1e-9)
        width = int(round(38 * frac))
        lines.append(f"  {label} {len(in_bin):>3} units  {'#' * width:<38}"
                     f" {100 * frac:>5.1f}% of perfusion")
    return "\n".join(lines)


def unit_detail(r: Result, index: int) -> str:
    units = list(r.units or [])
    if not units:
        return "  no per-unit data"
    if not 1 <= index <= len(units):
        return f"  there is no unit {index} (units run 1..{len(units)})"
    u = units[index - 1]
    co = max(r.co_l_min, 1e-9)
    q_tot = sum(x.perf for x in units) * co
    return (f"  unit {index}   z = {u.z:.2f} (0 non-dependent, 1 dependent)\n"
            f"    perfusion {u.perf * co:.2f} of {q_tot:.2f} L/min"
            f" ({100 * u.perf / max(sum(x.perf for x in units), 1e-9):.1f}% of flow)"
            f"   ventilation {u.va:.0f} mL/min   V/Q {u.vq:.3f}\n"
            f"    capillary PcCO2 {u.pc_co2:.1f} mmHg   PcO2 {u.po2:.1f} mmHg"
            f"   content {u.cc_o2:.3f} mL O2/mL blood\n"
            f"    disease {100 * u.diseased:.0f}%   open fraction"
            f" {100 * u.recruited:.0f}%   overdistension {100 * u.overdist:.0f}%"
            f"   {'unventilated (shunt)' if u.shunt else 'vented'}")


# --------------------------------------------------------------- session

@dataclass
class Session:
    scenario: Scenario | None = None
    sim: Simulation | None = None
    history: list[Result] = field(default_factory=list)

    def start(self, key: str | None = None) -> None:
        scen = BY_KEY[key] if key else self.scenario
        self.scenario = scen
        self.sim = scen.build() if scen is not None else Simulation(Patient(), Vent())
        self.history = [self.sim.current(keep_units=False)]

    def reset(self) -> None:
        key = self.scenario.key if self.scenario else None
        self.scenario = None
        self.start(key)

    # ------------------------------------------------------------ setters
    def apply(self, items: list[Assignment]) -> list[str]:
        assert self.sim is not None
        changed: list[str] = []
        for a in items:
            namespace, attr = a.key.split(".", 1)
            if namespace == "vent":
                self.sim.set_vent(**{attr: a.value})
            else:
                self.sim.set_patient(**{attr: a.value})
                if attr == "hco3":
                    # the metabolic set point: move the state with it
                    self.sim.hco3 = float(a.value)      # type: ignore[arg-type]
            changed.append(f"{attr}={a.value}")
        return changed

    # ------------------------------------------------------------ queries
    def settle_minutes(self, co_l_min: float) -> float:
        """Time for ~95% of the way to the operating point (tissue store)."""
        return 3.0 * 30.0 / max(co_l_min, 0.5)

    def abg_text(self) -> str:
        cur = self.sim.current()
        out = [f"  now       t = {cur.t_min:.1f} min"]
        out.append(gas_block(cur))
        out.append(flags_block(cur))
        note = effort_note(cur)
        if note:
            out.append(note)
        out.append(f"  settings  {settings_line(self.sim.vent, self.sim.patient)}")
        out.append(f"  patient   {patient_line(self.sim.patient)}")
        tgt = self.sim.target()
        out.append(f"  heading   pH {tgt.ph:.2f}   PaCO2 {tgt.paco2:.1f}"
                   f"   PaO2 {tgt.pao2:.1f}   HCO3 {tgt.hco3:.1f}"
                   f"   (about {self.settle_minutes(tgt.co_l_min):.0f} min away)")
        return "\n".join(out)

    def json_dump(self, history: list[Result] | None = None) -> str:
        cur = self.sim.current(keep_units=True)
        rows = history if history is not None else self.history

        def plain(r: Result) -> dict:
            d = asdict(r)
            d.pop("units", None)
            return d

        payload = {
            "scenario": self.scenario.key if self.scenario else None,
            "t_min": self.sim.t_min,
            "vent": asdict(self.sim.vent),
            "patient": {k: v for k, v in asdict(self.sim.patient).items()
                        if k != "lung"},
            "result": plain(cur),
            "units": [asdict(u) for u in (cur.units or [])],
            "history": [plain(r) for r in rows],
        }
        return json.dumps(payload, indent=1)

    def optimize_report(self, grid_tokens: list[str], top: int = 6,
                        weights: Weights | None = None) -> str:
        w = weights or Weights()
        p, v = self.sim.patient, self.sim.vent
        grids = grid_from_tokens(grid_tokens) if grid_tokens else None
        cur = self.sim.target()
        base_cost, base_parts = cost_of(p, v, cur, v, w)
        lines = [f"  current   {settings_line(v, p)}"
                 f"   pH {cur.ph:.2f} CO2 {cur.paco2:>4.0f} O2 {cur.pao2:>5.0f}"
                 f" Pplat {cur.pplat:>3.0f} dP {cur.dp_stat:>3.0f}"
                 f" MP {cur.mp_j_min:>3.0f} sh {cur.shunt_pct:>3.0f}%"
                 f"   cost {base_cost:.1f}  [{explain_from_parts(base_parts)}]"]
        cands = optimize(p, v, grids=grids, top=top, weights=w,
                        hco3=self.sim.hco3)
        lines.append(f"  ranked    {len(cands)} best of the search grid"
                     "  (lower cost = kinder to this lung)")
        for i, c in enumerate(cands, start=1):
            mark = "=" if same_settings(c.vent, v) else " "
            r = c.result
            lines.append(
                f"  {mark}{i}.  cost {c.cost:>7.1f}"
                f"   {settings_line(c.vent, p)}"
                f"   pH {r.ph:.2f} CO2 {r.paco2:>4.0f} O2 {r.pao2:>5.0f}"
                f" Pplat {r.pplat:>3.0f} dP {r.dp_stat:>3.0f}"
                f" MP {r.mp_j_min:>3.0f} sh {r.shunt_pct:>3.0f}%"
                f"   [{explain(c)}]")
        return "\n".join(lines)


def explain_from_parts(parts: dict, limit: int = 3) -> str:
    items = sorted(parts.items(), key=lambda kv: -kv[1])[:limit]
    return ", ".join(f"{k} {v:.0f}" for k, v in items) or "no penalties"


def same_settings(a: Vent, b: Vent) -> bool:
    keys = ["mode", "fio2", "peep", "vt_ml_kg", "rr", "ie", "dp"]
    return all(getattr(a, k) == getattr(b, k) for k in keys)


# --------------------------------------------------------------- shell

HELP = """\
  <blank>                 re-show the blood gas
  abg                     blood gas, mechanics, flags, where it is heading
  set k=v [k=v ...]       change ventilator or patient settings
  pset k=v ...            same, patient only
  vent | patient          print the current settings
  run <dur>               advance the clock (run 10m, run 2h)
  watch <dur> [every Nd]  print the time course as a table
  trend [n]               the last n recorded rows (default 20)
  vq                      per-unit ventilation/perfusion table
  vq <unit>               one unit in detail
  scenario [key]          describe or load a scenario
  scenarios               list the teaching scenarios
  optimize [grid]         rank settings, e.g. optimize peep=5,10,14 fio2=0.5,0.8
  reset                   back to the scenario's starting state
  help                    this text
  quit                    leave

  ventilator settings   mode  fio2  peep  vt (ml/kg of IBW)  rr  ie  dp|ps  ti
  patient settings      compliance  resistance  severity  scatter  hgb  co
                        sedation  muscle  renal  hco3  acid  bicarb_load  vco2
                        rq  temp  exp_limit  perf_gradient  vessels  pfo  hpv
                        n_units  height  weight  sex  name  notes
  durations             45s, 10m, 2h
"""


def do_command(session: Session, line: str) -> str | None:
    line = line.strip()
    if not line:
        return session.abg_text()
    low = line.lower()

    if low in ("q", "quit", "exit"):
        raise SystemExit(0)
    if low in ("h", "?") or low == "help" or low.startswith("help "):
        return HELP
    if low.startswith("scenarios") or low in ("scenario list", "scenario l"):
        return "\n".join(f"  {s.key:<20} {s.title}" for s in SCENARIOS)
    if low.startswith("scenario"):
        rest = line[len("scenario"):].strip()
        if not rest:
            s = session.scenario
            if s is None:
                return "  no scenario loaded; `scenario <key>` loads one"
            return (f"  {s.key}: {s.title}\n  {s.story}\n"
                    f"  try: {', '.join(s.try_)}\n  expect: {s.expect}")
        if rest.lower() not in BY_KEY:
            return f"  unknown scenario {rest!r}; `scenarios` lists them"
        session.start(rest.lower())
        return f"  loaded scenario {rest.lower()}\n" + session.abg_text()
    if low.startswith("reset"):
        session.reset()
        return "  reset to the starting state\n" + session.abg_text()
    if low == "vent" or low.startswith("vent "):
        return f"  {settings_line(session.sim.vent, session.sim.patient)}"
    if low == "patient" or low.startswith("patient "):
        return f"  {patient_line(session.sim.patient)}"
    if low.startswith("trend"):
        rest = line[len("trend"):].strip()
        n = int(rest) if rest.isdigit() else 20
        return trend_table(session.history[-n:])
    if low.startswith("vq"):
        rest = line[len("vq"):].strip()
        full = session.sim.current(keep_units=True)
        return (unit_detail(full, int(rest)) if rest.lstrip("-").isdigit()
                else vq_block(full))
    if low.startswith("abg"):
        return session.abg_text()
    if low.startswith("optimize"):
        rest = line[len("optimize"):].strip()
        return session.optimize_report(rest.split())
    if low.startswith("run"):
        rest = line[len("run"):].strip()
        minutes = parse_duration(rest or "10m")
        rows = time_to(session, minutes, every=max(minutes / 10.0, 0.1))
        final = rows[-1] if rows else session.sim.current()
        return (f"  t = {final.t_min:.1f} min\n" + gas_block(final)
                + "\n" + flags_block(final))
    if low.startswith("watch"):
        rest = line[len("watch"):].strip()
        every = None
        m = re.search(r"every\s+(\S+)", rest)
        if m:
            every = parse_duration(m.group(1))
            rest = rest[:m.start()].strip()
        minutes = parse_duration(rest or "1h")
        rows = time_to(session, minutes, every=every or max(minutes / 12.0, 0.5))
        return trend_table(rows)
    if low == "set" or low.startswith("set ") or low.startswith("pset"):
        patient_only = low.startswith("pset")
        rest = line[len("pset") if patient_only else len("set"):].strip()
        changed = session.apply(parse_assignments(rest.split(), patient_only=patient_only))
        return "  set " + ", ".join(changed) + "\n" + session.abg_text()
    if "=" in line:
        changed = session.apply(parse_assignments(line.split()))
        return "  set " + ", ".join(changed) + "\n" + session.abg_text()
    return f"  unknown command {line!r} (try `help`)"


def time_to(session: Session, minutes: float, every: float) -> list[Result]:
    """Run the clock, collecting a sample every `every` minutes."""
    dt = min(every / 4.0, 0.25)
    rows = session.sim.run(minutes, sample_every=every, dt=dt)
    session.history.extend(rows)
    return rows


def interactive(session: Session) -> int:
    try:
        import readline            # history on the arrow keys, if available
    except Exception:              # pragma: no cover
        pass
    print(f"ventsim {__version__} - ventilator and blood-gas simulation."
          "  `help` for commands, `quit` to leave.")
    print(session.abg_text())
    while True:
        try:
            line = input("\nvent> ")
        except EOFError:
            print()
            return 0
        except KeyboardInterrupt:
            print("  ^C - type quit to leave")
            continue
        try:
            out = do_command(session, line)
        except SystemExit:
            return 0
        except ValueError as exc:
            print(f"  {exc}")
            continue
        if out:
            print(out)


# --------------------------------------------------------------- entry point

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="ventsim",
        description="Multi-compartment ventilator and blood-gas simulation.")
    ap.add_argument("--scenario", metavar="KEY",
                    help="load a teaching scenario (see --list-scenarios)")
    ap.add_argument("--set", dest="set", nargs="*", default=[], metavar="K=V",
                    help="set ventilator/patient parameters")
    ap.add_argument("--run", metavar="DUR", help="advance the clock, then report")
    ap.add_argument("--watch", metavar="DUR", help="print a time table over DUR")
    ap.add_argument("--every", metavar="DUR", help="sampling interval for --watch")
    ap.add_argument("--optimize", nargs="*", metavar="GRID",
                    help="rank settings; GRID like peep=5,10,14 fio2=0.5,0.8")
    ap.add_argument("--script", metavar="FILE", help="run commands from a file")
    ap.add_argument("--json", action="store_true", help="machine-readable output")
    ap.add_argument("--list-scenarios", action="store_true")
    ap.add_argument("--quiet", action="store_true", help="suppress the ABG panel")
    args = ap.parse_args(argv)

    if args.list_scenarios:
        for s in SCENARIOS:
            print(f"{s.key:<20} {s.title}")
            print(f"{'':<20} {s.story}")
        return 0

    if args.scenario is not None and args.scenario not in BY_KEY:
        print(f"unknown scenario {args.scenario!r}; --list-scenarios lists them",
              file=sys.stderr)
        return 2

    session = Session()
    session.start(args.scenario)
    batch = bool(args.set or args.run or args.watch or args.json
                 or args.optimize is not None or args.script)

    try:
        if args.set:
            session.apply(parse_assignments(args.set))
        rows: list[Result] = []
        if args.watch:
            minutes = parse_duration(args.watch)
            every = parse_duration(args.every) if args.every else max(minutes / 12.0, 0.5)
            rows = time_to(session, minutes, every=every)
            print(trend_table(rows))
        elif args.run:
            minutes = parse_duration(args.run)
            rows = time_to(session, minutes, every=max(minutes / 10.0, 0.1))
            if not args.quiet:
                final = rows[-1]
                print(f"  t = {final.t_min:.1f} min")
                print(gas_block(final))
                print(flags_block(final))
        if args.optimize is not None:
            print(session.optimize_report(args.optimize))
        if args.json:
            print(session.json_dump(rows or None))
        if args.script:
            text = Path(args.script).read_text().splitlines()
            for ln in text:
                if not ln.strip() or ln.lstrip().startswith("#"):
                    continue
                out = do_command(session, ln)
                if out:
                    print(out)
        if not batch:
            return interactive(session)
    except SystemExit:
        return 0
    except (ValueError, OSError) as exc:
        print(f"  {exc}", file=sys.stderr)
        return 1
    if not args.quiet and not (rows or args.optimize is not None or args.json
                               or args.script):
        print(session.abg_text())
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
