# DEVELOPMENT — architecture, invariants, and how to work on this

This is the document to read before changing anything. The user-facing story
is in [README.md](README.md); where every physiological constant comes from is
in [REFERENCES.md](REFERENCES.md). This file is about the code: what each
piece does, which rules must not be broken, and how to prove you haven't
broken them.

**State as of last verified pass:** 74 pytest tests passing · 51/51 Python↔JS
parity cases at ~1e-14 · 12/12 scenario agreement · calibration clean ·
Playwright UI check worst delta 0.046, no JS errors. If any of that is not
true when you arrive, something regressed — the gates are in §4.

---

## 1. The shape of the code

```
                    ventsim/physiology.py      pure numerical primitives
                           │
        ventsim/lung.py ───┴─── ventsim/model.py
        (24 V/Q units)          (Vent, Patient, solve_steady → Result)
                           │
                    ventsim/dynamics.py         Simulation: state + time
                           │
     ventsim/scenarios.py ─┼─ ventsim/optimizer.py
     (12 teaching cases)   │    (penalty sum + coordinate descent)
                           │
                      ventsim/cli.py            interactive shell

     web/ventsim.js  ← a hand port of ALL of the above to JavaScript
     web/ventsim.html ← the UI, which loads only ventsim.js
```

### `physiology.py` (178 lines)
No state, no patient concepts. Henderson–Hasselbalch, the Siggaard-Andersen
base excess, the Hill O2 curve with Bohr/temperature shifts (`so2`, `p50`,
`o2_content`, `content_to_o2`), CO2 content (linear), the alveolar gas
equation, anion gap. Everything here is export-grade: both ports implement
every function.

### `lung.py` (344 lines)
`Lung` — one patient's anatomy as 24 compartments (apex→dependent `z`,
perfusion gradient, closing/overdistension pressure tables built once in
`build()`). `exchange(...)` — the gas exchange solve. CO2 in **closed form**:
`PaCO2 = S·K/(1−S)` with `S = Σ wᵢ/(1+Vᵢ/Qᵢ)`; this is exact because the CO2
content curve is linear. O2 by per-unit mass-balance bisection
(`_solve_unit_po2`, 26 steps) inside a damped mixed-venous-content outer loop
(12 passes, damping 0.55). Units with V/Q < 0.05 are identified as shunt once
per solve (`shunt_i[]`) and excluded from the vented loop.

### `model.py` (600 lines)
`Vent` (settings, `clamped()`), `Patient` (resolved defaults: IBW-scaled VCO2,
FRC, ve_rest), `ibw_kg` (ARDSNet), `drive_index`, `solve_steady` → `Result`
(~50 fields incl. flags as `(level, text)`), `result_flags`. The PSV
equilibrium lives here: `_drive_equilibrium` brackets `gap(d) = d −
drive_index(result gases)` then does exactly 20 bisection steps. A **pinned**
`drive=` argument gives the instantaneous operating point; `drive=None` gives
the self-consistent equilibrium.

### `dynamics.py` (207 lines)
`Simulation` — blood-gas state (`_fast/_slow` CO2 and O2 stores, `_hco3_base`
+ `_hco3_buf`, `drive`), advanced by **exact exponential relaxation**
(`x += (target−x)·(1−e^(−dt/τ))`), so any dt is stable *except* the PSV drive
loop, which `step()` sub-steps at `MAX_DRIVE_STEP_MIN = 0.1` min whenever
`mode == "psv" and sedation < 1` (§3, invariant 5). `adapt(hours)` runs 60 h
at dt 2.0 for renal compensation. `target(keep_units)` vs `current()`:
operating point vs where the blood actually is.

### `scenarios.py` (231 lines)
Twelve `Scenario` objects with `story` / `try_` / `expect` / `adapt_hours`.
`build()` returns a settled `Simulation`. The scenario *parameters* are
calibration, not measurement (see REFERENCES.md §"Authored, not cited").

### `optimizer.py` (274 lines)
`cost_of(patient, vent, result, base, weights)` → transparent penalty sum
(`hypoxaemia`, `hypercapnia`, `plateau`, `driving_pressure`, `local_stress`,
`fio2_dose`, `mechanical_power`, `work_of_breathing`, `support_deficit`,
`change`, …). `optimize(...)` — `method="grid"` (exhaustive, ≤60 000 combos)
or `"coord"` (default: two starts × single knobs + `JOINT_PAIRS`
`(vt_ml_kg,rr)`, `(peep,fio2)` × 8 sweeps; finds the cartesian optimum in
~1.4 s vs ~20 s). `parse_grid` parses `peep=5,8,12` token lists.

### `cli.py` (659 lines)
`parse_duration` / `parse_assignments` (alias maps for knobs) /
`grid_from_tokens`, renderers (`gas_block`, `vq_block`, `unit_detail`,
`trend_table`), the `Session` object, `do_command`, and `main()` with
`--scenario --set --run --watch --every --optimize --script --json
--list-scenarios --quiet`.

### `web/ventsim.js` (1 091 lines) and `web/ventsim.html` (961 lines)
The port must stay **numerically identical** to Python — that is the entire
point of the parity harness. The HTML builds controls dynamically into
`controls[group][field]`, has view buttons `data-v ∈ {steady, sim}` (labels
*steady* / *now*), and calls `sim.target(true)` / `sim.current(true)`.
Generated files `web/cases.json` and `web/scenarios.json` are **build
artifacts** — never hand-edit (§3, invariant 1).

### `tools/`
| tool | purpose | run when |
|---|---|---|
| `gen_cases.py` | regenerates `web/cases.json` (51 reference cases) | after **any** numeric change to `ventsim/` |
| `check_parity.js` | node: replays every case through `ventsim.js`, compares | after regenerating cases; also in pytest |
| `dump_scenarios.py` | regenerates `web/scenarios.json` | after changing `scenarios.py` |
| `check_scenarios.js` | node: the HTML's embedded SCENARIOS == Python's | after changing either side |
| `calibrate.py` | prints the tuning sweeps (no assertions) | sanity-checking scenario behaviour |
| `check_ui.py` | Playwright drives the real page, compares panel to Python | after touching `web/` at all; screenshots → `.artifacts/` |
| `probe_js.py` | dumps JS sim state for one scenario | debugging parity |

### `tests/` (7 files, 74 tests)
Behaviour tests, not frozen values: `test_physiology` (curve anchors),
`test_model` (mechanics/recruitment/PSV uniqueness/flags), `test_dynamics`
(time constants, renal, post-hypercapnic alkalosis, coarse-dt stability),
`test_scenarios` (clinical ranges), `test_optimizer` (parse, coord ≤ grid,
dedupe), `test_cli` (parsing + commands + `--json`), `test_parity` (runs the
node checks; skips without node).

---

## 2. Environment

* Python **3.14.4** (stdlib only — no third-party imports in `ventsim/`),
  Node **v22.22.1**, `pytest` and `playwright` installed, Chrome via
  `channel="chrome"` with `--no-sandbox --disable-gpu`.
* No `pyproject.toml` — the CLI is run from the repo root:
  `python3 -m ventsim …`.
* `gh` (snap) fails under the DSH sandbox with a DBus transient-scope error;
  it works with full process scope. Git push goes over HTTPS with gh's
  credential helper (`gh auth setup-git`, done once).
* `/tmp` is not reliable between tool calls — write artifacts into the
  workspace.
* Local preview server used in past sessions:
  `python3 -m http.server 8081 --bind 127.0.0.1 -d web` (8000/8080 were
  occupied on this machine).

---

## 3. Invariants — the rules that keep the two ports one simulator

1. **Every numeric change to `ventsim/*.py` must be mirrored in
   `web/ventsim.js`, then `web/cases.json` regenerated** (`python3
   tools/gen_cases.py web/cases.json`), then `node tools/check_parity.js`
   must exit 0. One-way doors: if you change Python and not JS, parity dies;
   if you regenerate cases from changed Python without changing JS, the same.
2. **Deterministic solves only.** O2 bisection 26 steps, drive bisection 20,
   mechanics 24 damped passes, content loop 12 passes at 0.55 — these exact
   counts are why JS reproduces Python to 1e-14. Never replace them with
   tolerance-based convergence or `scipy`.
3. **Golden-angle scatter sequences** (perfusion `mask`, `vent_het`) must be
   byte-identical across ports — same integer arithmetic, same order.
4. **Parity tolerances are not to be loosened**: |Δ| ≤ 0.02 abs or 2e-5 rel
   per field; `pao2/cao2/cvao2_diff/do2` ≤ 0.05; dynamic cases ≤ 0.05.
5. **The PSV sub-stepping** (`MAX_DRIVE_STEP_MIN`, mirrored in both ports)
   must survive any refactor. Without it, steps > ~0.5 min make the
   drive↔gases loop limit-cycle (PaCO2 swinging 22↔63 was the observed
   failure). The two coarse-step cases `course_psv_drive` (dt 5.0) and
   `course_psv_step_down` (dt 1.0) exist in cases.json to catch exactly this.
6. **HTML doctype stays first** — anything before it flips the page into
   quirks mode. (This is why per-file licence headers, if ever added, go
   *after* the doctype in `ventsim.html`.)
7. **Scenario baselines are load-bearing teaching numbers.** The table in
   README.md is generated from the live model; if a change moves a case out
   of its clinical range, the tests (`test_scenarios.py`) should fail —
   believe them, then re-tune the *patient*, not the test.

---

## 4. The verification gate (run before declaring anything done)

```bash
python3 -m pytest -q                                # 74 passed
node tools/check_parity.js                          # 51 cases, failing: 0, PARITY OK
python3 tools/dump_scenarios.py web/scenarios.json  # only if scenarios.py changed
node tools/check_scenarios.js web/scenarios.json    # scenarios agree (12 cases)
python3 tools/calibrate.py                          # rc 0, eyeball the sweeps
python3 tools/check_ui.py                           # needs playwright+chrome; delta < 0.05
```

Performance expectations on this machine: `solve_steady` VC ≈ 1.5 ms, PSV ≈
35 ms (drive bisection), optimizer coordinate descent ≈ 1.4 s, exhaustive
grid of 14 400 ≈ 20 s, `gen_cases.py` ≈ 30 s (it includes the PSV dynamic
courses), calibration ≈ 10 s, parity ≈ 2 s. A pytest run ≈ 38 s. If any of
these suddenly doubles, suspect an accidental O(n²) or a lost cache.

---

## 5. The git / GitHub state

* Remote: `https://github.com/tareko/ventilator-sim` (public, branch `main`,
  license AGPL-3.0 auto-detected). Identity: Tarek Loubani <tarek@tarek.org>.
* Commits: `714f4a8` initial upload · `348db05` REFERENCES.md.
* `.gitignore` keeps out `__pycache__`, `.pytest_cache`, `.artifacts/`.

## 6. Deliberately not done (the pick-up list)

* **`pip install -e .`** — no `pyproject.toml`; adding one plus a console
  script entry point is ~15 lines and untested.
* **Per-file SPDX headers** (`AGPL-3.0-or-later`) — offered, not yet done;
  touches all 30 files and needs the full gate re-run afterwards (and mind
  invariant 6 for the HTML).
* **A version tag** (`v0.1.0`) once you consider it pinned.
* Model ideas that were discussed but never built: oesophageal-pressure /
  transpulmonary-pressure display, a proper work-of-breathing Campbell
  diagram, breath-to-breath variability, paediatric patients (the IBW and
  FRC scalings are adult), V'CO2 changing with fever/agitation.
* `tools/check_ui.py` covers four scenarios (`ards_mod`, `weaning`,
  `chronic_hypercapnic`, `asthma`) — extending it is cheap and worthwhile if
  the UI grows new views.

## 7. Gotchas learned the hard way (all fixed; here so they stay fixed)

* `so2()` **does** apply the Bohr shift via `p50(ph, pco2, temp)` when no
  explicit P50 is passed — it silently ignored those arguments once, and the
  content bookkeeping disagreed with the unit solve.
* The pCO2 coefficient in `p50()` is **0.0010** per mmHg (carbamino effect
  only, ~0.6 mmHg P50 per 10 mmHg) — an earlier 0.0047 double-counted the
  Bohr effect and put P50 at 41 mmHg for a PaCO2 of 80.
* `anion_gap`'s albumin correction is `+2.5·(4.0 − albumin)` — a low albumin
  *hides* a gap.
* PSV semantics: `target()` pins the current drive (the machine's
  instantaneous answer); the equilibrium with a *free* drive is
  `solve_steady(..., drive=None)`. The dynamic run converges to the latter
  and they agree to 0.1 mmHg — a property the tests rely on.
* The optimizer needs the joint `(vt_ml_kg, rr)` move: the Vt↔RR cost
  surface is non-separable and single-knob descent parks at cost 119 where
  the joint move reaches 5.
* Shell-quoted inline `node -e` scripts break on quotes; put node checks in
  real files under `tools/` (that's why they exist as files).
