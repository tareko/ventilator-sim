"""Setting search: try a grid of settings and rank them by clinical cost.

The score is deliberately transparent - a weighted sum of interpretable
penalties rather than a magic number, so you can argue with it:

  hypoxaemia     SaO2 below the target band (or PaO2 below 62 mmHg)
  hyperoxia      PaO2 above 120 mmHg on enriched oxygen
  pH band        outside 7.30-7.50; worse the further out
  stress         driving pressure over 15, plateau over 28, peak local stress
  power          tidal mechanical power over ~12 J/min
  trapping       auto-PEEP over 3 cmH2O
  flow           cardiac output below the untreated baseline
  effort         neural drive above ~1.35x resting, and (in pressure support)
                 support that cannot deliver the tidal volume being demanded
  dose           a preference for low FiO2 once the target is met, and for
                 small changes from the settings you are already using

`python -m ventsim --scenario ards_mod --optimize` runs it with the defaults.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from itertools import product

from .model import Patient, Result, Vent, solve_steady


@dataclass
class Weights:
    hypoxaemia: float = 900.0
    hyperoxia: float = 6.0
    ph_low: float = 260.0
    ph_high: float = 220.0
    driving: float = 26.0
    plateau: float = 22.0
    stress: float = 8.0
    power: float = 3.0
    trapping: float = 9.0
    flow: float = 24.0
    rate: float = 2.0
    drive: float = 30.0
    fio2_dose: float = 1.2
    change: float = 3.0
    sao2_target: float = 0.92
    pao2_floor: float = 62.0
    pao2_ceiling: float = 120.0
    ph_lo: float = 7.30
    ph_hi: float = 7.50
    plateau_limit: float = 28.0
    driving_limit: float = 15.0
    stress_limit: float = 28.0
    power_limit: float = 12.0
    trapping_limit: float = 3.0
    rate_limit: float = 34.0
    drive_limit: float = 1.35


DEFAULT_GRIDS_VC = {
    "peep": [5, 8, 10, 12, 14, 16, 18, 20, 22, 24],
    "fio2": [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    "vt_ml_kg": [4.0, 5.0, 6.0, 7.0, 8.0],
    "rr": [12, 14, 16, 18, 20, 22, 24, 26, 28],
    "ie": [1.5, 2.0, 3.0, 4.0],
}
DEFAULT_GRIDS_PSV = {
    "peep": [5, 8, 10, 12, 14, 16, 18, 20],
    "fio2": [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    "dp": [6, 8, 10, 12, 14, 16, 18, 20, 24],
}

OPTIMISABLE = ("peep", "fio2", "vt_ml_kg", "rr", "dp", "ie")
# knob pairs that have to move together to be found at all
JOINT_PAIRS = [("vt_ml_kg", "rr"), ("peep", "fio2")]
# how far a knob has to move to count as one step of change
CHANGE_SCALE = {"peep": 3.0, "fio2": 0.1, "vt_ml_kg": 1.0, "rr": 3.0,
                "dp": 3.0, "ie": 1.0}


@dataclass
class Candidate:
    vent: Vent
    result: Result
    cost: float = 0.0
    parts: dict = field(default_factory=dict)

    def describe(self) -> str:
        v, r = self.vent, self.result
        if v.mode == "vc":
            mid = f"Vt {v.vt_ml_kg:.1f} RR {v.rr:.0f} I:E 1:{v.ie:g}"
        else:
            tag = "PS" if v.mode == "psv" else "dP"
            mid = f"{tag} {v.dp:.0f} RR {v.rr:.0f}"
        return (f"{v.label:<6} FiO2 {v.fio2:.2f} PEEP {v.peep:>2.0f}  {mid:<22}"
                f"| pH {r.ph:5.2f} PaCO2 {r.paco2:4.1f} PaO2 {r.pao2:5.1f} "
                f"Pplat {r.pplat:4.1f} dP {r.dp_stat:4.1f} aPEEP {r.auto_peep:4.1f} "
                f"CO {r.co_l_min:3.1f}  cost {self.cost:7.1f}")


def cost_of(patient: Patient, vent: Vent, res: Result, base: Vent,
            w: Weights, prefer_no_change: bool = True) -> tuple[float, dict]:
    p: dict[str, float] = {}

    def put(name: str, value: float) -> None:
        if value > 1e-9:
            p[name] = value

    low_sat = max(w.sao2_target - res.sao2, 0.0) / 0.02
    low_pa = max(w.pao2_floor - res.pao2, 0.0) / 8.0
    put("hypoxaemia", w.hypoxaemia * max(low_sat, low_pa))
    if res.fio2 > 0.30:
        put("hyperoxia", w.hyperoxia * max(res.pao2 - w.pao2_ceiling, 0.0))
    put("acidemia", w.ph_low * max(w.ph_lo - res.ph, 0.0) / 0.01)
    put("alkalemia", w.ph_high * max(res.ph - w.ph_hi, 0.0) / 0.01)
    put("driving_pressure", w.driving * max(res.dp_stat - w.driving_limit, 0.0))
    put("plateau", w.plateau * max(res.pplat - w.plateau_limit, 0.0))
    put("local_stress", w.stress * max(res.local_stress - w.stress_limit, 0.0))
    put("mechanical_power", w.power * max(res.mp_j_min - w.power_limit, 0.0))
    put("auto_peep", w.trapping * max(res.auto_peep - w.trapping_limit, 0.0))
    put("cardiac_output", w.flow * max(0.80 * patient.co_l_min - res.co_l_min, 0.0) / 0.1)
    put("high_rate", w.rate * max(res.rr - w.rate_limit, 0.0))
    # Work of breathing: only bites when the patient is awake enough to have a
    # neural drive, which is exactly when under-supporting them matters.
    put("work_of_breathing", w.drive * max(res.drive - w.drive_limit, 0.0) / 0.25)
    if vent.mode == "psv":
        put("support_deficit", w.drive * 0.6 * max(1.0 - res.effort_ratio, 0.0) / 0.05)
    if low_sat <= 0:
        put("fio2_dose", w.fio2_dose * max(res.fio2 - 0.30, 0.0) * 10.0)
    if prefer_no_change:
        delta = 0.0
        for key, scale in CHANGE_SCALE.items():
            a, b = getattr(base, key), getattr(vent, key)
            if a is not None and b is not None and scale > 0:
                delta += abs(b - a) / scale
        put("change", w.change * delta)
    return sum(p.values()), p


def optimize(
    patient: Patient,
    vent: Vent | None = None,
    *,
    grids: dict | None = None,
    knobs: list[str] | None = None,
    weights: Weights | None = None,
    top: int = 6,
    prefer_no_change: bool = True,
    hco3: float | None = None,
    method: str = "auto",
) -> list[Candidate]:
    """Rank settings for this patient.

    With an explicit `grids` (see `parse_grid`) every combination is tried.
    Otherwise coordinate descent walks the default grids one knob at a time:
    the same neighbourhood in about a hundredth of the solves, which matters
    when this is sitting behind a keypress.
    """
    vent = (vent or Vent()).clamped()
    w = weights or Weights()
    explicit = grids is not None
    if grids is None:
        grids = dict(DEFAULT_GRIDS_PSV if vent.mode == "psv" else DEFAULT_GRIDS_VC)
    if knobs is None:
        knobs = [k for k in grids if grids.get(k)]
    else:
        knobs = [k for k in knobs if grids.get(k)]
    for k in knobs:
        if k not in OPTIMISABLE:
            raise ValueError(f"cannot optimise {k!r}; choose from {OPTIMISABLE}")
    if not knobs:
        raise ValueError("nothing to search: no knobs with a grid")

    def evaluate(v: Vent) -> Candidate:
        res = solve_steady(patient, v, hco3=hco3) if hco3 is not None \
            else solve_steady(patient, v)
        cost, parts = cost_of(patient, v, res, vent, w, prefer_no_change)
        return Candidate(v, res, cost, parts)

    if method == "auto":
        method = "grid" if explicit else "coord"

    if method == "grid":
        combos: list[dict] = [{}]
        for k in knobs:
            combos = [dict(c, **{k: val}) for c in combos for val in grids[k]]
        if len(combos) > 60000:
            raise ValueError(f"grid too large ({len(combos)} combinations);"
                             " ask for fewer values or use method='coord'")
        out: list[Candidate] = []
        for combo in combos:
            try:
                v = replace(vent, **combo).clamped()
            except (ValueError, TypeError):
                continue
            out.append(evaluate(v))
    else:
        # Coordinate descent over the default grids.  Single knobs alone are
        # not enough here: dropping Vt only pays if the rate rises with it, so
        # each sweep also tries the pairs that interact, from two starts.
        moves = [(k,) for k in sorted(knobs)]
        moves += [(a, b) for a, b in JOINT_PAIRS if a in knobs and b in knobs]

        def descend(start: Vent) -> list[Candidate]:
            local: list[Candidate] = []
            best = evaluate(start)
            local.append(best)
            current = best.vent
            for _sweep in range(8):
                improved = False
                for names in moves:
                    grids_for = [grids[n] for n in names]
                    for tup in product(*grids_for):
                        kw = {n: t for n, t in zip(names, tup)}
                        if all(getattr(current, n) == t for n, t in kw.items()):
                            continue
                        try:
                            cand = evaluate(replace(current, **kw).clamped())
                        except (ValueError, TypeError):
                            continue
                        local.append(cand)
                        if cand.cost < best.cost - 1e-9:
                            best, current, improved = cand, cand.vent, True
                if not improved:
                    break
            return local

        if vent.mode == "psv":
            alt = replace(vent, dp=14.0, rr=14.0, peep=max(vent.peep, 12.0))
        else:
            alt = replace(vent, vt_ml_kg=6.0, rr=20.0, peep=max(vent.peep, 12.0))
        pool: list[Candidate] = []
        for start in (vent, alt.clamped()):
            pool.extend(descend(start))
        out = pool
    unique: dict[tuple, Candidate] = {}
    for c in out:
        key = signature(c.vent)
        prev = unique.get(key)
        if prev is None or c.cost < prev.cost:
            unique[key] = c
    out = list(unique.values())
    out.sort(key=lambda c: c.cost)
    return out[:max(top, 1)]


def signature(vent: Vent) -> tuple:
    return tuple(getattr(vent, k) for k in
                 ("mode", "fio2", "peep", "vt_ml_kg", "rr", "ie", "dp"))


def parse_grid(spec: str) -> dict:
    """'peep=5,8,12 fio2=0.4,0.6' -> {'peep': [5,8,12], 'fio2': [0.4,0.6]}"""
    grids: dict[str, list[float]] = {}
    for chunk in str(spec).split():
        key, sep, vals = chunk.partition("=")
        if not sep:
            raise ValueError(f"expected key=v1,v2 - got {chunk!r}")
        try:
            nums = [float(x) for x in vals.split(",") if x.strip()]
        except ValueError:
            raise ValueError(f"bad value list in {chunk!r}") from None
        if not nums:
            raise ValueError(f"no values for {key!r}")
        grids.setdefault(key.strip().lower(), []).extend(nums)
    if not grids:
        raise ValueError("empty grid specification")
    return grids


def explain(c: Candidate, limit: int = 4) -> str:
    if not c.parts:
        return "no penalties"
    items = sorted(c.parts.items(), key=lambda kv: -kv[1])[:limit]
    return ", ".join(f"{k} {v:.0f}" for k, v in items) or "no penalties"
