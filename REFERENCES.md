# References — where the physiology in this simulator comes from

Every number in the model is either taken from the literature, a textbook
standard, or was **calibrated by us** to reproduce the textbook behaviour of a
disease. That last class matters: a teaching simulator is not a clinical trial,
and pretending every constant has a citation would be dishonest. So this file
labels every entry:

* **[L]** — *literature-pinned.* A formula or constant with a specific source
  (section [Bibliography](#bibliography)).
* **[T]** — *textbook-standard.* The kind of number found in any respiratory
  physiology text (West, Levitsky, Nunn); cited at book level rather than to a
  specific paper.
* **[C]** — *calibrated/authored.* Chosen by us, tuned in
  `tools/calibrate.py` until the scenarios behaved like their names. These are
  the numbers to challenge first if a case ever stops teaching the right
  lesson.

Code references are to `ventsim/physiology.py` unless another file is named.

---

## 1. Acid–base chemistry

| quantity | value | source | |
|---|---|---|---|
| Henderson–Hasselbalch equation | pH = 6.101 + log10(HCO3 / (0.0301·PCO2)) | pKa for the CO2/HCO3 system at 37 °C, and whole-blood CO2 solubility α = 0.0301 mmol·L⁻¹·mmHg⁻¹; standard values in any acid–base text (Siggaard-Andersen 1974, §Bibliography) | **[L]** |
| base excess | BE = (1 − 0.014·Hb)·[HCO3 − 24 + (1.43·Hb + 7.7)·(pH − 7.4)] | The Van Slyke equation as published by **Siggaard-Andersen** (*The Acid-Base Status of the Blood*, 4th ed. 1974, p 51). Worked examples were cross-checked against three independent calculators during development: [pediatriconcall.com](https://www.pediatriconcall.com/calculators/base-excess-calculator), the [base-excess review in PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC1297616/), and labtestsguide.com | **[L]** |
| anion gap, albumin-corrected | (Na − Cl − HCO3) + 2.5·(4.0 − albumin) | the standard "add 2.5 per g/dL of albumin below 4" correction **[T]** | **[T]** |

## 2. Oxygen carriage

| quantity | value | source | |
|---|---|---|---|
| O2–Hb curve | Hill equation with n = 2.7, P50 = 26.6 mmHg | Hill (1910) formulated the equation precisely for haemoglobin's O2 curve (see [Wikipedia: Hill equation](https://en.wikipedia.org/wiki/Hill_equation_(biochemistry))); n ≈ 2.7 and P50 ≈ 26.6 mmHg for adult HbA are the standard figures (e.g. [ScienceDirect overview](https://www.sciencedirect.com/topics/medicine-and-dentistry/oxygen-haemoglobin-dissociation-curve), citing Severinghaus 1966) | **[L]** |
| curve shifts | log10 P50 = 0.44·(7.4 − pH) + 0.0010·(PCO2 − 40) + 0.024·(T − 37) | Bohr effect coefficient −0.44 pH-units from the *in vitro* whole-blood literature summarised by **Severinghaus's 1966 nomogram**; the small independent CO2 term (~0.6 mmHg of P50 per 10 mmHg PCO2) is the carbamino effect *after* the pH shift is counted; ~5.7 % per °C for temperature | **[L]** |
| Hb O2 capacity | 1.34 mL O2 per g Hb (Hüfner's number) | **[T]** | **[T]** |
| dissolved O2 | 0.00003 mL O2 per mL blood per mmHg (Henry's law) | **[T]** | **[T]** |
| CO2 content slope | 0.0121 mmol·L⁻¹·mmHg⁻¹, taken as linear over the physiological range | a straight-line fit to the standard whole-blood CO2 dissociation curve; the linearity is what makes the closed-form multi-compartment CO2 solution possible (§4) | **[C]** |

## 3. Gas pressures and the alveolar gas equation

| quantity | value | source | |
|---|---|---|---|
| barometric pressure, water vapour | PB 760 mmHg, PH2O 47 mmHg at 37 °C | **[T]** | **[T]** |
| PiO2 | FiO2 × (PB − PH2O) | **[T]** | **[T]** |
| alveolar gas equation | PAO2 = PiO2 − PaCO2·(FiO2 + (1−FiO2)/RQ), RQ = 0.8 | the full form that reduces correctly to PiO2 − PaCO2/RQ on room air; **[T]** (West, *Respiratory Physiology*, gas exchange chapter) | **[T]** |
| BTPS→STPD conversion | (713/760)·(273/15/310.15) ≈ 0.826 | **[T]** | **[T]** |

## 4. The lung: compartments, V/Q, recruitment, shunt

The architecture — a vertical stack of compartments, perfusion rising towards
the dependent zones, ventilation rising less, so V/Q falls from apex to base —
is **West's**, and the multi-compartment V/Q description of gas exchange with
its low-V/Q, high-V/Q and shunt modes follows his Chapter 5 (*Respiratory
Physiology: The Essentials*, 9th–11th ed.; also Levitsky, *Pulmonary
Physiology*, Ch. 5). The model's 24 units and their distributions are our
discretisation of it.

| quantity | value | source | |
|---|---|---|---|
| gravity perfusion gradient | linear apex→dependent, ratio set per patient | West's zone model, **[T]** | **[T]** |
| anatomical + bronchial-venous shunt | 2 % of cardiac output in a normal lung | the classic ~2 % "anatomical shunt" figure **[T]** | **[T]** |
| anatomical dead space | 2.2 mL per kg predicted body weight | **[T]** (≈150 mL in a 70 kg adult) | **[T]** |
| FRC | 35 mL/kg IBW | **[T]** (≈2.4 L at 70 kg) | **[T]** |
| CO2 in closed form | PaCO2 = S·K/(1−S), S = Σ wᵢ/(1+Vᵢ/Qᵢ) over compartments, K = V'CO2/(β·Q̇) | our derivation, using the linear CO2 content curve; it is the standard multi-compartment mass balance done exactly instead of iteratively | **[C]** |
| Vd/Vt (Enghoff) | (PaCO2 − mixed expired PCO2)/PaCO2 | **Enghoff (1938)** replaced alveolar PCO2 with arterial PCO2 in Bohr's equation; in general use since (see [dead-space review, PMC4857382](https://pmc.ncbi.nlm.nih.gov/articles/PMC4857382/)) | **[L]** |
| hypoxic pulmonary vasoconstriction | perfusion redirected away from low-V/Q units with a saturating gain | **von Euler & Liljestrand (1946)** — the observation; our gain curve is calibrated (see [Swenson's review](https://en.wikipedia.org/wiki/Hypoxic_pulmonary_vasoconstriction)) | concept **[L]**, curve **[C]** |
| closing pressures of diseased units | recruit between 7 and 30 cmH2O | the ARDS recruiting range observed in PEEP-titration studies; the specific spread is our calibration | **[C]** |
| overdistension stress limits | ~17 cmH2O onset in healthy units, higher in remodelled ones | the notion that regional transpulmonary pressure, not airway pressure, injures is **Gattinoni's** stress/strain programme; the numeric limits are calibrated | concept **[L]**, numbers **[C]** |
| ventilation heterogeneity | scatter of ventilation between units, widened by disease | **[T]** concept, **[C]** magnitude | **[C]** |

## 5. Mechanics and ventilator-induced harm

| quantity | value | source | |
|---|---|---|---|
| driving pressure | Pplat − total PEEP, flagged above 15, hard-warning above 17 cmH2O | **Amato et al., NEJM 2015** — driving pressure was the ventilation variable that best stratified survival in ARDS ([PubMed 25693014](https://pubmed.ncbi.nlm.nih.gov/25693014/)) | **[L]** |
| tidal mechanical power | 0.098·RR·[½·E·Vt² + R·Vt²/Ti] J/min | **Gattinoni et al., Intensive Care Med 2016** ("Ventilator-related causes of lung injury: the mechanical power"); the ≥17 J/min attention threshold follows that programme | **[L]** |
| auto-PEEP | trapped volume above FRC against the expiratory time constant, τ = R·C | the *existence and detection* of intrinsic PEEP in obstructive disease is **Pepe & Marini, Am Rev Respir Dis 1982** ([PubMed 7046541](https://pubmed.ncbi.nlm.nih.gov/7046541/)); our one-compartment calculation is textbook **[T]** | **[L]** |
| P/F ratio and its strata | PaO2/FiO2; mild/moderate/severe ARDS at 201–300 / 101–200 / ≤100 | **Berlin definition, JAMA 2012** ([PubMed 22797452](https://pubmed.ncbi.nlm.nih.gov/22797452/)) | **[L]** |
| oxygenation index | OI = 100·FiO2·Pmean/PaO2 | standard paediatric critical-care metric, as used in the **Pediatric Acute Lung Injury Consensus (2015)** definition ([PubMed 25647235](https://pubmed.ncbi.nlm.nih.gov/25647235/)) | **[L]** |
| predicted body weight | male 50 + 0.91·(h − 152.4); female 45.5 + … | exactly as published in the **ARDSNet ARMA trial, NEJM 2000** ([full text](https://www.nejm.org/doi/full/10.1056/NEJM200005043421801)) | **[L]** |
| lung-protective defaults (6 mL/kg, Pplat ≤ 30) | | the same ARDSNet trial | **[L]** |

## 6. Haemodynamics

Cardiac output falls with mean airway pressure and with overdistension
(capillary compression). The *phenomenon* is textbook — intrathoracic pressure
impeding venous return, PEEP reducing output in ARDS — but our specific
slopes (output starting to fall above a mean pressure of 14 cmH2O, losing
~2.5 % per further cmH2O, more with `rv_sensitivity` raised for a struggling
right ventricle) are **[C]**: calibrated so that "PEEP 22 halves the cardiac
output" comes out as it does at the bedside. There is no Starling curve here.

## 7. Time: stores, buffering, renal compensation

| quantity | value | source | |
|---|---|---|---|
| two gas stores per gas | lung+blood τ ≈ 0.6·FRC/V̇A (CO2); tissue τ = 30/CO min (CO2), 10/CO (O2) | the *multi-compartment store* picture (a fast well-mixed compartment and a slow tissue compartment) is the classic CO2-stores description of Farhi and colleagues; the specific time constants are our simplification, chosen so blood gases settle in minutes and tissue equilibration takes ~half an hour at resting flow | concept **[T]**, constants **[C]** |
| acute metabolic compensation | ~+0.10 mmol/L HCO3 per mmHg PaCO2, τ = 36 s | the classic acute respiratory-acidosis rule "bicarbonate rises ~1 mmol per 10 mmHg"; buffering itself is near-instant, we give it half a minute for the tissues | rule **[T]**, τ **[C]** |
| chronic renal compensation | ~+0.35 mmol/L per mmHg (≈3.5 per 10), τ = 12 h, and *asymmetric* — losing bicarbonate uses a 0.10/mmHg slope | chronic respiratory acidosis expects ~+3.5–4 mmol/L per 10 mmHg, with the renal response underway by 6–12 h and maximal over days (see e.g. [anaesthesiamcq.com, §4.5](https://www.anaesthesiamcq.com/AcidBaseBook/ab4_5.php)). One 12-hour time constant is a deliberate simplification, and the asymmetry is our encoding of post-hypercapnic alkalosis being slow to resolve | rule **[T]**, τ and asymmetry **[C]** |
| bicarbonate distribution space | 0.5 L/kg | the classic ~40 % body-water figure used when correcting acidosis | **[T]** |

## 8. Spontaneous breathing and respiratory drive

The controller is a classic chemostat: a PaCO2 threshold (set point ~34 mmHg),
a linear hypercapnic ventilatory response, a hypoxic term that only matters
when saturation falls, both blunted by sedation, with a neural lag. The *shape*
is textbook (West; Ganong; the Duffin/Whitelaw control literature), but our
slopes are deliberately tame:

| quantity | value | note | |
|---|---|---|---|
| hypercapnic response | 1.05 L/min per mmHg above threshold | the awake normal response is usually quoted steeper (2–3 L/min/mmHg); ours is a sedated-ICU-patient response | **[C]** |
| hypoxic response | +30 % V'E at SaO2 0.82 | **[C]** | **[C]** |
| neural lag | τ = 0.4 min | central chemoreceptor dynamics act over tens of seconds; we round up | **[C]** |
| pressure support equilibrium | bracketed bisection on drive | the fixed-point iteration oscillates (loop gain ≈ −22), which is why we solve it as a root-finding problem and sub-step the time course at 0.1 min — numerical method, not physiology | **[C]** |

## 9. Reference ranges and flag thresholds

The ABG panel's warning thresholds — pH 7.25/7.55, PaO2 60, PaCO2 50, Pplat
30, driving pressure 15–17, mechanical power 17 J/min, shunt 20 %, local
stress 35 cmH2O — are the standard rule-of-thumb targets of protective
ventilation as taught from the ARDSNet/Amato/Gattinoni programme above. They
are judgement calls, and `ventsim/model.py:result_flags` is where to argue
with them.

---

## Bibliography

**Papers**

1. Amato MBP, Meade MO, Slutsky AS, et al. *Driving pressure and survival in
   the acute respiratory distress syndrome.* **N Engl J Med** 2015;372:747–55.
   [PubMed 25693014](https://pubmed.ncbi.nlm.nih.gov/25693014/)
2. Gattinoni L, Tonetti T, Cressoni M, et al. *Ventilator-related causes of
   lung injury: the mechanical power.* **Intensive Care Med** 2016;42:1567–75.
   (Companion review: Gattinoni L, et al. *Mechanical power and potential for
   lung protection in ARDS.* Intensive Care Med 2018;44:856–63.)
3. Pepe PE, Marini JJ. *Occult positive end-expiratory pressure in
   mechanically ventilated patients with airflow obstruction: the auto-PEEP
   effect.* **Am Rev Respir Dis** 1982;126:166–70.
   [PubMed 7046541](https://pubmed.ncbi.nlm.nih.gov/7046541/)
4. The ARDS Network (Brower RG, Matthay MA, et al.). *Ventilation with lower
   tidal volumes as compared with traditional tidal volumes for acute lung
   injury and the acute respiratory distress syndrome.* **N Engl J Med**
   2000;342:1301–8. [Full text](https://www.nejm.org/doi/full/10.1056/NEJM200005043421801)
   — the source of the predicted-body-weight formula and the 6 mL/kg / Pplat-30
   defaults.
5. ARDS Definition Task Force (Ranieri VM, et al.). *Acute respiratory
   distress syndrome: the Berlin definition.* **JAMA** 2012;307:2526–33.
   [PubMed 22797452](https://pubmed.ncbi.nlm.nih.gov/22797452/)
6. Enghoff H. *Volumen inefficax: Bemerkungen zur Frage des schädlichen
   Raumes.* **Uppsala Läkarefören Förh** 1938;44:191–218 — the Enghoff
   modification of Bohr's dead space. (Context:
   [dead-space assessment in ARDS, PMC4857382](https://pmc.ncbi.nlm.nih.gov/articles/PMC4857382/))
7. von Euler US, Liljestrand G. *Observations on the pulmonary arterial blood
   pressure in the cat.* **Acta Physiol Scand** 1946;12:301–20 — hypoxic
   pulmonary vasoconstriction. (Review:
   [Swenson, High Alt Med Biol 2013](https://en.wikipedia.org/wiki/Hypoxic_pulmonary_vasoconstriction))
8. Siggaard-Andersen O. *The Acid-Base Status of the Blood*, 4th revised ed.
   Munksgaard, Copenhagen, 1974 — the base-excess equation (p 51). (Worked
   examples cross-checked against
   [pediatriconcall.com](https://www.pediatriconcall.com/calculators/base-excess-calculator),
   [PMC1297616](https://pmc.ncbi.nlm.nih.gov/articles/PMC1297616/), and
   labtestsguide.com during calibration.)
9. Hill AV. *The possible effects of the aggregation of the molecules of
   haemoglobin on its dissociation curves.* **J Physiol** 1910;40(Suppl):i–vii
   — the Hill equation, written for this exact curve. (Context:
   [Wikipedia: Hill equation](https://en.wikipedia.org/wiki/Hill_equation_(biochemistry)))
10. Severinghaus JW. *Blood gas calculator.* **J Appl Physiol** 1966;21:1108–16
    — the standard nomogram for the Bohr/temperature shifts of the O2
    dissociation curve.
11. Pediatric Acute Lung Injury Consensus Conference Group (Khemani RG, Smith
    LS, Zimmerman JJ, et al.). *Pediatric acute respiratory distress syndrome:
    definition, incidence, and epidemiology.* **Pediatr Crit Care Med**
    2015;16(Suppl 1):S17–34. [PubMed 25647235](https://pubmed.ncbi.nlm.nih.gov/25647235/)
    — the oxygenation index as a severity metric.
12. Kallet RH, et al., on the mis-estimates inherent in dead-space formulas —
    why Enghoff Vd/Vt (which we use) reads high when shunt is present.
    [Intensive Care Med 2013](https://www.sciencedirect.com/science/article/abs/pii/S1569904813002267)

**Books**

* West JB. *Respiratory Physiology: The Essentials.* 9th–11th ed. Wolters
  Kluwer. — V/Q distributions, zones, gas exchange; the intellectual skeleton
  of `ventsim/lung.py`.
* Levitsky MG. *Pulmonary Physiology.* 9th ed. McGraw-Hill, 2018 — the same
  material from another angle (Ch. 5, ventilation–perfusion relationships).
* Nunn's *Applied Respiratory Physiology.* 8th ed. Elsevier — time constants,
  dead space, closing volume, the details of the mechanics layer.
* Siggaard-Andersen's monograph above, for everything base excess.
* Lumb AB. *Nunn's Applied Respiratory Physiology.* (same book; listed for the
  acid-base chapters as well.)

**Web resources used during calibration** (worked examples, not authorities):

* [anaesthesiamcq.com, Acid-Base §4.5](https://www.anaesthesiamcq.com/AcidBaseBook/ab4_5.php)
  — renal compensation rules for chronic respiratory acidosis.
* [LITFL: Driving pressure](https://litfl.com/driving-pressure/) and
  [LITFL: ARDSNet ventilation strategy](https://litfl.com/ardsnet-ventilation-strategy/)
  — quick statements of the Amato and ARDSNet positions.
* [Deranged Physiology: the p50 value of a blood gas sample](https://derangedphysiology.com/main/cicm-primary-exam/respiratory-system/Chapter-1131/p50-value-blood-gas-sample)
  — P50 = 26.6 mmHg and its use.

---

## Authored, not cited — the honest list

Nothing below has a source. Each was chosen because it makes the model behave
the way the textbooks say the disease behaves, and each is a simplification a
reviewer should feel free to challenge (they live in `ventsim/lung.py`,
`ventsim/model.py` and `ventsim/dynamics.py`):

* the 24-compartment discretisation itself, and the golden-angle scatter used
  to distribute ventilation heterogeneity between units;
* `affected_fraction = 0.55·severity` — how much of the lung a given
  "severity" takes out;
* closing-pressure and overdistension-limit spreads per unit
  (`7–30 cmH2O`, `17/30 cmH2O` onsets), and `stiff_affected = 0.7`;
* the tidal-swing share of recruitment (`RECRUIT_TIDAL = 0.35`);
* the right-ventricle/haemodynamic slopes of §6;
* the drive slopes and lags of §8;
* the tissue-store time constants (30/CO and 10/CO minutes);
* the scenario patients themselves — every compliance, resistance, severity
  and heterogeneity in `ventsim/scenarios.py` is a calibration aimed at a
  textbook description, tuned in `tools/calibrate.py`.

If you want to check any of them, that is exactly what `tools/calibrate.py`
is for: it prints the sweeps those numbers were tuned against.
