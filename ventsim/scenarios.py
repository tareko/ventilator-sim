"""Bedside teaching cases: who is in the bed, and what is worth trying.

Each scenario carries the story and the commands worth typing, so the tool can
teach itself.  `build()` returns a ready `Simulation`; scenarios marked with
`adapt_hours` have already had their kidneys catch up (chronic CO2 retention).
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .dynamics import Simulation
from .model import Patient, Vent


@dataclass
class Scenario:
    key: str
    title: str
    patient: Patient
    vent: Vent
    story: str = ""
    try_: list[str] = field(default_factory=list)
    expect: str = ""
    adapt_hours: float = 0.0

    def build(self) -> Simulation:
        """A simulation sitting at this scenario's settings, gases at baseline."""
        sim = Simulation(self.patient, self.vent)
        if self.adapt_hours > 0:
            sim.adapt(self.adapt_hours)
            sim.equilibrate()
        return sim


def _s(key, title, patient_kw, vent_kw, story, try_, expect, adapt_hours=0.0) -> Scenario:
    return Scenario(key, title, Patient(**patient_kw), Vent(**vent_kw),
                    story, list(try_), expect, adapt_hours)


SCENARIOS: list[Scenario] = [
    _s(
        "normal", "Healthy adult, anaesthetised",
        dict(name="healthy 42 y", height_cm=172, weight_kg=75, compliance=60,
             resistance=5, sedation=1.0, notes="reference case"),
        dict(mode="vc", vt_ml_kg=8, rr=12, peep=5, fio2=0.30),
        "Normal lungs under general anaesthesia. Use this as the ruler: every "
        "other case is judged by how far it is from these numbers.",
        ["abg", "set rr=20", "run 3m", "set rr=12 vt_ml_kg=12", "run 3m", "watch 5m"],
        "PaCO2 follows the ventilation change within a minute or two; the A-a "
        "gradient stays small; plateau pressure is ~14 cmH2O.",
    ),
    _s(
        "ards_mild", "Mild ARDS",
        dict(name="ARDS, mild", height_cm=170, weight_kg=70, severity=0.28,
             compliance=42, resistance=7, heterogeneity=0.40, hgb=12.0,
             co_l_min=5.5, rv_sensitivity=1.1,
             notes="P/F 200-280, both lungs hazy"),
        dict(mode="vc", vt_ml_kg=6, rr=19, peep=8, fio2=0.50),
        "Early ARDS: mostly non-dependent collapse, still a good deal of "
        "aerated lung. Protective ventilation should be comfortable here.",
        ["abg", "optimize", "set vt_ml_kg=8", "abg"],
        "Driving pressure stays under 15; oxygenation responds to modest PEEP. "
        "Going to 8 ml/kg looks attractive and is still wrong for the lung, "
        "not the numbers.",
    ),
    _s(
        "ards_mod", "Moderate ARDS",
        dict(name="ARDS, moderate", height_cm=170, weight_kg=70, severity=0.55,
             compliance=28, resistance=8, heterogeneity=0.45, hgb=11.0,
             co_l_min=6.0, rv_sensitivity=1.2,
             notes="P/F ~100, dependent collapse"),
        dict(mode="vc", vt_ml_kg=6, rr=21, peep=10, fio2=0.80),
        "The classic ARDS trade-off case: dependent lung is collapsed, the "
        "non-dependent lung is already near its limit. PEEP is a negotiation "
        "between oxygenation, stress and cardiac output.",
        ["abg", "optimize", "run 10m", "set peep=20", "run 10m", "trend"],
        "Recruitment climbs with PEEP while shunt falls, but local stress and "
        "dead space rise and cardiac output drops. Somewhere in the middle "
        "driving pressure is lowest - that is the honest optimum.",
    ),
    _s(
        "ards_severe", "Severe ARDS",
        dict(name="ARDS, severe", height_cm=168, weight_kg=68, severity=0.80,
             compliance=18, resistance=10, heterogeneity=0.55, hgb=9.5,
             co_l_min=7.0, rv_sensitivity=1.5,
             notes="P/F <100, tiny baby lung"),
        dict(mode="vc", vt_ml_kg=4.5, rr=32, peep=14, fio2=1.0),
        "Almost nothing left to ventilate. Every millilitre is expensive and "
        "every centimetre of PEEP costs cardiac output. This is where you "
        "start thinking about proning, paralysis and ECMO rather than knobs.",
        ["abg", "set vt_ml_kg=6", "abg", "optimize"],
        "6 ml/kg produces a plateau pressure over 40. CO2 cannot be cleared "
        "without destroying the lung: accept the hypercapnia.",
    ),
    _s(
        "pneumonia", "Lobar pneumonia",
        dict(name="right lower lobe pneumonia", height_cm=175, weight_kg=78,
             severity=0.30, affected_fraction=0.42, heterogeneity=0.15,
             compliance=45, resistance=6, hgb=12.0, co_l_min=6.5,
             lung_overrides=dict(affected_closing_lo=12.0, affected_closing_hi=26.0),
             notes="dense consolidation, little V/Q scatter"),
        dict(mode="vc", vt_ml_kg=6, rr=22, peep=10, fio2=0.70),
        "A big mass of lung that is perfused and not ventilated: pure shunt, "
        "plus a bit of recruitable collapse around it.",
        ["abg", "set fio2=1.0", "abg", "set peep=18", "run 5m", "abg"],
        "FiO2 to 1.0 moves PaO2 much less than you expect - that is the "
        "definition of shunt. PEEP helps some, at the price of plateau "
        "pressure and cardiac output.",
    ),
    _s(
        "copd", "COPD / emphysema on the ventilator",
        dict(name="COPD, emphysema", height_cm=170, weight_kg=70, severity=0.15,
             compliance=45, resistance=16, exp_limit=3.0, heterogeneity=0.45,
             vascular_loss=0.35, co_l_min=6.5, hgb=14.0, vco2=190.0, sedation=1.0,
             notes="obstruction + capillary loss"),
        dict(mode="vc", vt_ml_kg=8, rr=11, ie=2.0, peep=5, fio2=0.30),
        "Airflow obstruction with emphysematous capillary loss. Exhalation "
        "takes time, and time is what fast rates take away from it.",
        ["abg", "set rr=20", "run 3m", "abg", "set rr=10 ie=4.0", "run 3m", "trend"],
        "Faster rate: auto-PEEP climbs, plateau climbs, CO2 barely improves. "
        "Slower rate with a long exhalation: auto-PEEP falls and the same "
        "ventilation is achieved with less pressure.",
    ),
    _s(
        "asthma", "Severe asthma",
        dict(name="asthma, near-fatal", height_cm=178, weight_kg=75, severity=0.10,
             compliance=55, resistance=30, exp_limit=6.0, heterogeneity=0.50,
             co_l_min=7.0, vco2=210.0, sedation=1.0, notes="bronchospasm"),
        dict(mode="vc", vt_ml_kg=8, rr=12, ie=2.0, peep=5, fio2=0.35),
        "The same settings that are merely untidy in COPD kill in asthma: the "
        "expiratory time constant is enormous and breath stacking is quick.",
        ["abg", "set rr=8 ie=4.0", "run 5m", "trend", "set rr=16"],
        "At RR 12 the plateau is over 40 with auto-PEEP near 25 and cardiac "
        "output down a third. Slow the rate, take the I:E to 1:4 and accept "
        "the CO2 - permissive hypercapnia.",
    ),
    _s(
        "pe", "Massive pulmonary embolism",
        dict(name="PE, massive", height_cm=175, weight_kg=80, severity=0.05,
             compliance=55, resistance=6, heterogeneity=0.60,
             perfusion_gradient=0.8, vascular_loss=0.60, co_l_min=4.0,
             rv_sensitivity=2.0, vco2=230, sedation=1.0, pfo=0.06,
             notes="ventilated but not perfused"),
        dict(mode="vc", vt_ml_kg=8, rr=16, peep=5, fio2=0.40),
        "The lung is nearly normal and the vessels are not. Ventilation goes "
        "where blood does not, so nothing is exchanged and the expired CO2 "
        "looks tidy while the arterial CO2 does not.",
        ["abg", "set peep=15", "abg", "optimize"],
        "High Vd/Vt with nearly normal shunt, hypoxaemia that FiO2 partly "
        "fixes, and PEEP making the cardiac output worse - the dead-space "
        "pattern, opposite to ARDS.",
    ),
    _s(
        "fibrosis", "Idiopathic pulmonary fibrosis",
        dict(name="IPF", height_cm=165, weight_kg=60, severity=0.45,
             affected_fraction=0.95, compliance=22, resistance=5,
             heterogeneity=1.00, perfusion_gradient=1.6, co_l_min=6.0, hgb=13.0,
             lung_overrides=dict(healthy_closing=2.0, affected_closing_lo=2.0,
                                 affected_closing_hi=6.0, stiff_affected=0.40,
                                 od_lo_affected=26.0, od_span_affected=10.0),
             notes="stiff, not recruitable"),
        dict(mode="vc", vt_ml_kg=6, rr=30, peep=8, fio2=0.60),
        "Fibrotic lung is stiff but it is *open*: there is nothing to recruit "
        "and plenty to overdistend. The ARDS habit of pushing PEEP fails here.",
        ["abg", "set peep=16", "run 5m", "abg"],
        "PEEP barely moves PaO2 and moves dead space and stress up. Fast shallow "
        "breathing with low volumes suits this lung; the low compliance is not "
        "an invitation to bigger tidal volumes.",
    ),
    _s(
        "obesity", "Obese patient, post-operative",
        dict(name="obese 45 y", height_cm=170, weight_kg=145, compliance=30,
             resistance=8, frc_l=1.30, sedation=1.0, co_l_min=6.5, vco2=240,
             hgb=14.0,
             lung_overrides=dict(healthy_closing=18.0, od_lo_healthy=20.0),
             notes="heavy chest wall, low FRC, dependent atelectasis"),
        dict(mode="vc", vt_ml_kg=8, rr=14, peep=5, fio2=0.50),
        "A heavy chest wall and a low FRC: the lung closes off at ordinary "
        "filling pressures. This is the one case where high PEEP is the whole "
        "answer.",
        ["abg", "set peep=10", "run 3m", "set peep=14", "run 3m", "trend"],
        "PaO2 at PEEP 5 is bad because the dependent lung is simply collapsed; "
        "PEEP 10-14 recruits it and oxygenation jumps, with little cost because "
        "there is no inflammation to make the lung fragile.",
    ),
    _s(
        "chronic_hypercapnic", "Chronic CO2 retention (48 h)",
        dict(name="COPD, chronic hypercapnic", height_cm=168, weight_kg=65,
             severity=0.25, compliance=42, resistance=14, exp_limit=2.5,
             vco2=200, sedation=1.0, hgb=15.0, co_l_min=6.0,
             notes="kidneys fully compensated"),
        dict(mode="vc", vt_ml_kg=6, rr=13, ie=3.0, peep=6, fio2=0.35),
        "Two days of stable hypercapnia: the kidneys have grown the bicarbonate "
        "that makes the pH look acceptable. Now somebody decides to 'normalise' "
        "the PaCO2.",
        ["abg", "set rr=32", "run 30m", "abg", "set rr=13", "run 12h", "trend"],
        "pH goes from ~7.37 to ~7.68 in half an hour: post-hypercapnic "
        "alkalosis. Undoing it takes days - the kidney that grew the "
        "bicarbonate over two days does not forget in thirty minutes.",
        adapt_hours=48.0,
    ),
    _s(
        "weaning", "Spontaneous breathing, insufficient support",
        dict(name="day 4, weaning trial", height_cm=172, weight_kg=72,
             severity=0.30, compliance=45, resistance=7, sedation=0.15,
             muscle=0.8, vco2=240, hgb=11.0, co_l_min=5.5,
             notes="breathing spontaneously on PSV"),
        dict(mode="psv", dp=8, rr=8, peep=5, fio2=0.35),
        "Awake enough to drive, supported too little to satisfy the drive. "
        "Pressure support does not fix minute ventilation on its own; the "
        "patient does, and here the patient has to overwork.",
        ["abg", "set dp=16", "abg", "run 10m", "set dp=6", "run 10m", "trend"],
        "More support lowers the drive and the work; less support raises both. "
        "Unlike controlled modes, there is a real equilibrium per support level "
        "and the blood gas moves to it over minutes.",
    ),
]

BY_KEY = {s.key: s for s in SCENARIOS}


def get(key: str) -> Scenario:
    if key not in BY_KEY:
        raise KeyError(f"unknown scenario {key!r}; try one of: "
                       + ", ".join(BY_KEY))
    return BY_KEY[key]


def keys() -> list[str]:
    return list(BY_KEY)
