# ventsim — a bedside ventilator & blood-gas simulator

Turn the knobs on a ventilator and watch the arterial blood gases answer back —
including **how long the answer takes**. That last part is the point: a PEEP
change shows up in the PaO2 within minutes, a respiratory-rate change in the
PaCO2 in minutes, and a bicarbonate change over days.

Two ways to use it:

* a **command-line session** (`python3 -m ventsim`) with an ABG panel, a V/Q
  histogram, a trend table, a time clock and a setting search;
* a **single-file web page** (`web/ventsim.html`) — same numbers, dependency-free
  JavaScript port of the model, no server, no build step.

Nothing here is connected to a real ventilator and nothing here should be used
to make a clinical decision. It is a teaching instrument: the values are the
textbook ones, and the interesting part is *why* they move.

---

## Quick start

```bash
python3 -m ventsim --scenario ards_mod          # interactive session
python3 -m ventsim --list-scenarios             # what the cases teach
python3 -m ventsim --scenario copd --set rr=20 ie=3 --run 30m
python3 -m ventsim --scenario ards_mod --set peep=16 fio2=0.6 --optimize
python3 -m ventsim --scenario chronic_hypercapnic --watch 12h --every 2h
python3 -m ventsim --scenario pe --run 20m --json > pe.json
```

Open `web/ventsim.html` in a browser (double-click it, `file://` is fine) or
serve the folder: `python3 -m http.server -d web 8000`.

In the interactive session, `help` lists the commands. The useful ones:

| command | what it does |
|---|---|
| `peep=14 fio2=0.7 vt=6 rate=24 ps=12 ie=3` | change the ventilator, show the new operating point immediately |
| `pset severity=0.9 hgb=8 co=4 sedation=0.3` | change the *patient* instead |
| `run 10m` / `watch 2h every 15m` | let the clock go; `watch` prints a row per interval |
| `abg` `trend 8` `vq` `vq 7` `flags` | the panel, the last 8 samples, the 24 compartments, one compartment, just the warnings |
| `optimize` / `optimize vt=4,5,6 rr=16,20,24 peep=5,10,15` | rank settings by a transparent penalty sum |
| `scenario pneumonia` `reset` `settle 5m` | load a case, undo everything, wait for the gases to catch up |
| `script file.txt` | run commands from a file |

Durations are `45s`, `10m`, `2h`. Anything you set is clamped to the range a
real machine would let you into (`vt_ml_kg` 3–15, `rr` 4–60, `peep` 0–40,
`fio2` 0.21–1.0).

---

## What the model is

Four layers, each one a plain physiological argument, no hidden magic.

### 1. The lung: 24 compartments with a V/Q distribution

Gas exchange happens in 24 units spread from apex to dependent back, with a
West-style gravity gradient on perfusion and a distribution of closing and
overdistension pressures. A diseased fraction (severity-driven) gets low
compliance, higher closing pressure and a lower overdistension limit.

* **Recruitment.** A unit opens when the pressure above its closing pressure
  lasts long enough; PEEP therefore lowers shunt — up to the point where it
  stops helping and starts costing.
* **Overdistension.** Each unit has a local stress limit; the worst unit's
  stress is reported separately from the average, because in a heterogeneous
  lung it is the worst units that get hurt. `local stress` routinely sits 10
  cmH2O above the driving pressure.
* **Hypoxic pulmonary vasoconstriction** dials perfusion down towards
  badly-ventilated units; vascular loss (emboli, emphysema) removes capillary
  bed and shows up as dead space.
* **CO2 is solved in closed form.** With a linear CO2 content curve, mixed
  venous CO2 content and the arterial tension can be solved exactly for a
  multi-compartment lung: `PaCO2 = S·K/(1−S)`, where `S` sums the `1/(1+V/Q)`
  terms. Dead space (Enghoff), shunt and the A-a gradient all fall out of it.
* **O2 is solved per unit** by mass balance against the haemoglobin curve
  (Hill, n = 2.7, P50 26.6 mmHg with Bohr and temperature shifts), then mixed
  with the shunted blood. Unperfused units hold inspired gas; unventilated ones
  pass blood through untouched.

### 2. Mechanics

Compliance and resistance give you PIP, plateau, driving pressure, mean airway
pressure, peak flow and — where exhalation is incomplete — **auto-PEEP**. The
asthma case generates 25 cmH2O of it at a rate of 12 with I:E 1:2, which is
why its plateau pressure looks so frightening and its cardiac output so poor.
Tidal mechanical power is reported in J/min
(`0.098 · RR · [½·E·Vt² + R·Vt²/Ti]`).

### 3. Blood and buffers

Arterial content comes from the Hb curve; oxygen content, A-a gradient, shunt
fraction, OI, P/F, oxygen delivery and the Hb-corrected base excess are all
computed, never looked up. Base excess is the Van Slyke / Siggaard-Andersen
form:

```
BE = (1 − 0.014·Hb) · [HCO3 − 24 + (1.43·Hb + 7.7)·(pH − 7.4)]
```

Bicarbonate has two mechanisms with two different clocks:

* **acute cellular buffering**, τ = 36 s, about +0.1 mmol/L per mmHg of PaCO2;
* **renal compensation**, τ = 12 h (scaled by the patient's `renal` function),
  about +0.35 mmol/L per mmHg of chronic hypercapnia — and deliberately
  *asymmetric*: getting rid of bicarbonate again is slower than retaining it,
  which is what produces post-hypercapnic alkalemia when you ventilate a
  chronically retaining patient too hard.

### 4. Time

PaCO2 and PaO2 relax through **two stores**: a fast lung+blood store
(τ ≈ 0.6·FRC/V̇A for CO2) and a slow tissue store (τ = 30/CO minutes for CO2,
10/CO for O2). So the answer to a settings change arrives in two phases:
a quick move in the first minute or two, then a slow tail as the tissues
catch up. Every relaxation is exact (`x += (target − x)·(1 − e^{−dt/τ})`), so
the simulation is stable at any step size.

### 5. A patient who breathes back: pressure support

In `psv` mode the patient decides. Respiratory drive follows the classic
response line (`set_point` 34 mmHg, hypoxic gain, blunted by `sedation`),
relaxes with a 24 s neural lag, and drives tidal volume through the muscles
versus the pressure support the machine provides. The equilibrium is found by
bracketed bisection — the fixed-point iteration oscillates (loop gain ≈ −22) —
and the time course is sub-stepped to one tenth of the neural lag so the
drive↔gases loop can't limit-cycle. The panel reports `drive` (× resting
minute ventilation) and `effort_ratio` (how much of the work the support is
actually doing for them).

Two things worth knowing about those:

* The operating point the panel prints first **holds the current drive fixed**.
  That is the machine's instantaneous answer. To see the patient respond to a
  change in support, let the clock run (`run 30m`) — then the drive moves to
  where the blood gases put it, and the two agree exactly with the equilibrium
  the bisection finds. In the web page the same distinction is the difference
  between the *steady* and *now* views.
* On the `weaning` case: at PS 4 the patient is working hard and still
  failing — drive 1.6× resting, effort 66%, PaCO2 47. At PS 12 the effort is
  entirely gone (effort 1.00) and PaCO2 42. **Above PS 12 nothing changes at
  all** — no further fall in PaCO2, no further rise in tidal volume — however
  much support you add, because the patient has already taken the wheel.

---

---

## The teaching cases

```bash
python3 -m ventsim --list-scenarios      # full story, experiments and expectations
```

Steady state at each case's own settings:

| case | what it is | vent settings | pH | PaCO2 | PaO2 | HCO3 | Pplat | dP | aPEEP | Vd/Vt | shunt | CO | P/F |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `normal` | Healthy adult, anaesthetised | VC 8 ml/kg · RR 12 · PEEP 5 · FiO2 0.30 · I:E 1:2 | 7.38 | 41 | 135 | 24.0 | 14 | 9 | 0.0 | 0.24 | 2% | 5.0 | 450 |
| `ards_mild` | Mild ARDS | VC 6 ml/kg · RR 19 · PEEP 8 · FiO2 0.50 · I:E 1:2 | 7.39 | 41 | 117 | 24.0 | 19 | 11 | 0.0 | 0.34 | 14% | 5.5 | 233 |
| `ards_mod` | Moderate ARDS | VC 6 ml/kg · RR 21 · PEEP 10 · FiO2 0.80 · I:E 1:2 | 7.42 | 38 | 122 | 24.0 | 29 | 19 | 0.0 | 0.37 | 21% | 5.4 | 152 |
| `ards_severe` | Severe ARDS | VC 4 ml/kg · RR 32 · PEEP 14 · FiO2 1.00 · I:E 1:2 | 7.36 | 43 | 163 | 24.0 | 39 | 25 | 0.0 | 0.51 | 18% | 4.2 | 163 |
| `pneumonia` | Lobar pneumonia | VC 6 ml/kg · RR 22 · PEEP 10 · FiO2 0.70 · I:E 1:2 | 7.44 | 37 | 58 | 24.0 | 24 | 14 | 0.0 | 0.37 | 37% | 6.4 | 83 |
| `copd` | COPD / emphysema on the ventilator | VC 8 ml/kg · RR 11 · PEEP 5 · FiO2 0.30 · I:E 1:2 | 7.40 | 40 | 112 | 24.0 | 20 | 12 | 2.5 | 0.25 | 7% | 6.3 | 372 |
| `asthma` | Severe asthma | VC 8 ml/kg · RR 12 · PEEP 5 · FiO2 0.35 · I:E 1:2 | 7.42 | 39 | 141 | 24.0 | 45 | 16 | 24.5 | 0.29 | 2% | 2.8 | 403 |
| `pe` | Massive pulmonary embolism | VC 8 ml/kg · RR 16 · PEEP 5 · FiO2 0.40 · I:E 1:2 | 7.38 | 42 | 115 | 24.0 | 15 | 10 | 0.0 | 0.45 | 7% | 4.0 | 288 |
| `fibrosis` | Idiopathic pulmonary fibrosis | VC 6 ml/kg · RR 30 · PEEP 8 · FiO2 0.60 · I:E 1:2 | 7.54 | 29 | 125 | 24.0 | 27 | 19 | 0.0 | 0.41 | 16% | 5.8 | 208 |
| `obesity` | Obese patient, post-operative | VC 8 ml/kg · RR 14 · PEEP 5 · FiO2 0.50 · I:E 1:2 | 7.39 | 41 | 78 | 24.0 | 37 | 32 | 0.0 | 0.28 | 23% | 6.1 | 157 |
| `chronic_hypercapnic` | Chronic CO2 retention (48 h) | VC 6 ml/kg · RR 13 · PEEP 6 · FiO2 0.35 · I:E 1:3 | 7.37 | 55 | 119 | 30.7 | 17 | 10 | 0.8 | 0.33 | 10% | 6.0 | 340 |
| `weaning` | Spontaneous breathing, insufficient support | PSV PS 8 · RR 8 · PEEP 5 · FiO2 0.35 · I:E 1:2 | 7.36 | 44 | 85 | 24.0 | 15 | 10 | 0.0 | 0.33 | 14% | 5.5 | 244 |

Things worth an hour with each case:

* **`normal` is the ruler.** Note the A-a gradient of 28 at FiO2 0.30 and a
  shunt of 2%. Everything else is judged by how far it is from these numbers.
* **`ards_mod`: find the PEEP that helps.** PaO2 climbs as shunt closes, the
  plateau climbs with it, and cardiac output falls. `optimize` picks the
  trade-off and tells you which penalties it weighed.
* **`pneumonia` is a shunt, not a ventilation problem.** PaO2 60 at FiO2 0.70
  with a normal pH and PaCO2, shunt 37%. Push FiO2 to 1.0 and PaO2 moves to 70
  while the A-a gradient explodes from 400 to 600 — the signature of blood
  that never sees gas. PEEP is the only thing that helps: at PEEP 18 the shunt
  is down to 16% and PaO2 is 122, paid for with a plateau of 31 cmH2O.
* **`asthma` is air trapping.** Look at auto-PEEP of 24.5 and a cardiac output
  of 2.8. Slow it down and give it time to exhale (`rr=6 ie=4`): auto-PEEP
  falls to 8.4 and the output recovers to 5.4, Pplat 45 becomes 25 — and the
  PaCO2 climbs to 72 because minute ventilation was halved. That trade-off is
  the whole art of the case.
* **`pe` is dead space, not shunt.** Vd/Vt 0.45 with a shunt of 7%: the lung is
  ventilating normally into alveoli the blood isn't reaching. PEEP makes it
  worse in every direction — at PEEP 15 the Vd/Vt is up to 0.47, PaO2 has
  fallen from 115 to 87 and the cardiac output from 4.0 to 2.8.
* **`chronic_hypercapnic` has been like this for two days.** HCO3 30.7 and BE
  +4.6 at pH 7.37. Now set `rr=32` and `run 20m`: the PaCO2 falls, the
  bicarbonate does not, and the pH lands in the 7.7s. That is the reason you
  do not normalise a chronic retainer's PaCO2.
* **`weaning` is a patient on support that is too small.** Set `ps=4` and
  `run 30m`: drive 1.6× resting, effort 66%, and a PaCO2 of 47 despite the
  effort. Set `ps=12` and the work is gone. Any support beyond that changes
  nothing at all — same tidal volume, same gases — because the patient has
  already taken the wheel.

---

## Where the numbers come from

Every formula in the model is a standard one. Values quoted as used here.

| thing | formula / value | where it is from |
|---|---|---|
| pH, HCO3 | Henderson–Hasselbalch, pKa 6.101, CO2 solubility 0.0301 mmol/L/mmHg | standard; any acid-base text |
| base excess | Van Slyke / Siggaard-Andersen, Hb-corrected (see above) | Siggaard-Andersen's original blood-buffer line; worked examples cross-checked against [pediatriconcall.com](https://www.pediatriconcall.com/calculators/base-excess-calculator) and the [base-excess review in PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC1297616/) |
| anion gap | Na − Cl − HCO3, corrected 2.5 × (4.0 − albumin) | standard albumin correction |
| O2–Hb curve | Hill equation, n = 2.7, P50 = 26.6 mmHg | [Wikipedia: oxygen–hemoglobin dissociation curve](https://en.wikipedia.org/wiki/Oxygen%E2%80%93hemoglobin_dissociation_curve); n ≈ 2.7 is the usual fit |
| curve shifts | log10 P50 = 0.44·(7.4−pH) + 0.0010·(PCO2−40) + 0.024·(T−37) | *in vitro* Bohr coefficient for whole blood; the CO2 term is the carbamino effect *after* the pH effect is counted (~0.6 mmHg of P50 per 10 mmHg PCO2); ~5.7 % per °C |
| O2 content | 1.34 mL O2/g Hb × Hb × SaO2 + 0.00003 × PO2 (mL O2 per mL blood) | Hüfner number + dissolved term |
| CO2 content | 0.0121 mmol/L per mmHg, linear in tension | used for the closed-form CO2 solution |
| inspired O2 | PiO2 = FiO2 × (PB − 47) | water vapour at 37 °C |
| alveolar gas | PAO2 = PiO2 − PaCO2·(FiO2 + (1−FiO2)/RQ), RQ = 0.8 | alveolar gas equation |
| Vd/Vt | Enghoff modification: (PaCO2 − mixed-expired CO2)/PaCO2 | so anatomical dead space is included |
| OI, P/F | 100 × FiO2 × mean airway pressure / PaO2; PaO2/FiO2 | Berlin-definition arithmetic |
| ideal body weight | 50 (male) / 45.5 (female) + 0.91 × (height − 152.4) | ARDSNet formula |
| gravity gradient | perfusion scaled apex→dependent, V/Q falling with height | West's zone/V/Q distributions |
| driving pressure | Pplat − PEEP, target < 15 cmH2O; stress limit per compartment | Gattinoni's driving-pressure and stress/strain work |
| mechanical power | 0.098 · RR · [½·E·Vt² + R·Vt²/Ti], J/min | Gattinoni's tidal mechanical power |
| auto-PEEP | from incomplete exhalation: the trapped volume against the expiratory time constant (τ = R·C) | standard mechanics |
| FRC | 35 mL/kg IBW, anatomical dead space 2.2 mL/kg IBW | standard per-weight volumes |
| stores | lung+blood τ ≈ 0.6·FRC/V̇A (CO2), tissue τ = 30/CO min (CO2), 10/CO min (O2) | the multi-compartment argument for why blood gases take minutes and bicarbonate takes days |
| drive | response line about PaCO2 (set point 34 mmHg) plus a hypoxic term, blunted by sedation | classic CO2 response-line physiology |

The scenario parameters themselves (compliance, resistance, severity,
heterogeneity, vascular loss, PFO) were **tuned against the textbook
descriptions** of each disease: the numbers in `ventsim/scenarios.py` are the
patient, and the outputs in the table above come from running them.
`tools/calibrate.py` prints the sweeps used while choosing those parameters —
PEEP curves, tidal volume sweeps, Vd/Vt and shunt across the range of cases —
and is the file to look at if a case ever stops behaving like its name.

---

## How it is tested

```bash
python3 -m pytest -q                 # 74 tests: physiology, model, time course, CLI, parity
node tools/check_parity.js           # every Python number vs the JavaScript port
node tools/check_scenarios.js        # the web page's scenarios are the same patients
python3 tools/calibrate.py           # the tuning sweeps (no assertions, prints tables)
python3 tools/check_ui.py            # Playwright: drives the real web UI, compares numbers
```

* `tests/test_physiology.py` pins the acid-base and dissociation numbers.
* `tests/test_model.py`, `test_dynamics.py`, `test_scenarios.py` test the
  behaviour rather than exact values: recruitment goes the right way, PEEP
  eventually costs cardiac output, a coarse time step does not destabilise the
  drive loop, a chronic retainer arrives compensated, over-support stops
  mattering.
* `tests/test_parity.py` runs the JavaScript port against 51 Python-generated
  reference cases (`web/cases.json`, regenerated by `tools/gen_cases.py`).
  Tolerances are tight on purpose: |Δ| ≤ 0.02 absolute or 2e-5 relative per
  field (0.05 for the oxygen content family). In practice the two ports agree
  to 1e-14, including the time-stepped spontaneous-breathing courses.
* `tools/check_ui.py` needs `playwright` and Chrome; it skips nothing, it just
  drives the page the way you would and compares the printed panel to Python.

## Limitations, honestly

* **A teaching model, not a patient simulator.** Gas exchange is linear-ish,
  the CO2 content curve is a straight line, chest-wall and lung elastance are
  lumped into one compliance, and there is one heart with one number.
* Haemodynamics are a rule, not a solution: cardiac output falls with mean
  airway pressure scaled by `rv_sensitivity`. It reproduces "a PEEP of 22
  halves your cardiac output" but it is not a Starling curve.
* Recruitment is a pressure-time rule, not surfactant physics. HPV is a
  ventilatory-to-perfusion gain, not a vascular-smooth-muscle model.
* No alveolar ventilation/perfusion mismatch from the airways in spontaneous
  breathing beyond resistance, no breath-to-breath variability, no cough, no
  secretion, no asynchrony.
* The renal compensation term is one time constant plus an asymmetry, so it
  gets chronic compensation approximately right and never gives you a proper
  mixed-disorder puzzle unless you load acids and bicarbonate explicitly.

## Layout

```
ventsim/          physiology.py lung.py model.py dynamics.py scenarios.py optimizer.py cli.py
web/              ventsim.html  ventsim.js  cases.json  scenarios.json
tools/            gen_cases.py check_parity.js check_scenarios.js calibrate.py check_ui.py ...
tests/            pytest suite
```

## License

Every file in this repository is licensed under the **GNU Affero General
Public License v3.0** (AGPL-3.0). The full text is in [LICENSE](LICENSE), and
also at [gnu.org](https://www.gnu.org/licenses/agpl-3.0.html).

AGPL rather than GPL on purpose: it is the version of the copyleft that
applies when the software runs as a network service. If you run a modified
copy of this simulator somewhere people use it through the network, the source
of that modified copy has to be available to them.

Copyright &copy; 2026 Tarek Loubani.
