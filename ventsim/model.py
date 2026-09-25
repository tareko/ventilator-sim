"""Ventilator + patient model: mechanics, recruitment, dynamics, blood gases.

One operating point (`solve_steady`)
-----------------------------------
1.  Breath timing.  Volume control fixes tidal volume; pressure control fixes
    the insufflation pressure.  Pressure support lets the patient set the
    pattern (see `drive_index`).
2.  Recruitment and overdistension.  Every unit opens when the distending
    pressure beats its own closing pressure and is overdistended when local
    stress passes its own limit.  Both feed back on compliance, overdistension
    also strangles local perfusion, so this runs as a damped iteration.
3.  Ventilation is distributed by local compliance; perfusion follows gravity
    and is pulled *away* from unventilated units by hypoxic vasoconstriction.
4.  The V/Q solver (lung.py) returns PaO2, PaCO2, shunt, dead space.
5.  Indices: plateau/driving pressure, peak local stress, tidal mechanical
    power, OI, P/F, Enghoff dead-space fraction, base excess.

Time stepping (`ventsim.dynamics.Simulation`)
---------------------------------------------
Blood gases relax to the operating point through two compartments (fast
lung/blood store, slow tissue store); bicarbonate moves with acute cellular
buffering in minutes, renal compensation over hours, plus an acid or bicarbonate
load.  That is what makes the classic bedside stories appear on their own.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

from . import physiology as phys
from .lung import Lung, exchange

MODES = ("vc", "pc", "psv")

ANATOMIC_VD_PER_KG = 2.2      # mL dead space per kg IBW (154 mL at 70 kg)
FRC_PER_KG = 0.035            # L per kg IBW  -> ~2.5 L at 70 kg
OD_WIDTH = 12.0               # stress range over which a unit becomes overdistended
OD_STIFFEN = 0.5              # fraction of compliance lost at maximal overdistension
OD_CAPILLARY = 0.75           # fraction of perfusion lost at maximal overdistension
NONDEP_TRANSMISSION = 0.5     # PEEP is felt more by non-dependent units
OD_CO = 0.30                  # cardiac output lost per unit of mean overdistension
OD_MEAN_BASE = 14.0           # mean airway pressure above which CO starts to fall
OD_MEAN_CO = 0.25             # ... and this much per 10 cmH2O above it
RECRUIT_TIDAL = 0.35          # tidal swing's share in keeping units open, besides PEEP
DRIVE_SLOPE = 1.05            # L/min extra minute ventilation per mmHg above threshold
HYPOXIC_GAIN = 0.30           # hypoxic drive: fraction of V'E added at SaO2 0.82
PMUS_MAX = 12.0               # cmH2O extra muscle pressure per resting-drive unit
PMUS_CAP = 18.0               # ceiling on the muscle contribution
ACUTE_HCO3_PER_MMHG = 0.10    # cellular buffering, mmol/L per mmHg PaCO2
CHRONIC_HCO3_PER_MMHG = 0.35  # renal compensation, mmol/L per mmHg PaCO2
CHRONIC_HYPO_PER_MMHG = 0.10  # ... limited when PaCO2 is below normal
TAU_ACUTE_MIN = 0.6           # minutes
TAU_RENAL_MIN = 720.0         # minutes (12 h to reach most of the renal response)
BUFFER_SPACE_FRACTION = 0.5   # litres of distribution per kg body weight


def ibw_kg(height_cm: float, sex: str = "m") -> float:
    """Ideal body weight, ARDSNet formula (kg)."""
    if str(sex).lower().startswith("f"):
        return 45.5 + 0.91 * (height_cm - 152.4)
    return 50.0 + 0.91 * (height_cm - 152.4)


def _clamp(x: float, lo: float, hi: float) -> float:
    return lo if x < lo else hi if x > hi else x


# ---------------------------------------------------------------------------
# settings & patient
# ---------------------------------------------------------------------------
@dataclass
class Vent:
    """Ventilator settings.  `ie` is the expiratory-to-inspiratory ratio (2.0 = 1:2)."""

    mode: str = "vc"
    fio2: float = 0.5
    peep: float = 5.0
    vt_ml_kg: float = 8.0        # volume-control target, mL/kg IBW
    rr: float = 12.0             # set rate; in PSV this is the apnoea backup rate
    ie: float = 2.0
    dp: float = 12.0             # insufflation pressure above baseline (PC / PSV), cmH2O
    ti_s: float | None = None    # explicit inspiratory time, overrides ie

    def clamped(self) -> "Vent":
        v = replace(self)
        v.mode = str(v.mode).lower().strip()
        if v.mode not in MODES:
            raise ValueError(f"mode must be one of {MODES}, got {self.mode!r}")
        v.fio2 = _clamp(v.fio2, 0.21, 1.0)
        v.peep = _clamp(v.peep, 0.0, 45.0)
        v.vt_ml_kg = _clamp(v.vt_ml_kg, 2.0, 15.0)
        v.rr = _clamp(v.rr, 1.0, 60.0)
        v.ie = _clamp(v.ie, 0.35, 10.0)
        v.dp = _clamp(v.dp, 0.0, 60.0)
        if v.ti_s is not None:
            v.ti_s = _clamp(v.ti_s, 0.3, 4.0)
        return v

    @property
    def label(self) -> str:
        return {"vc": "VC-AC", "pc": "PC-AC", "psv": "PSV"}[self.mode]

    def summary(self, ibw: float) -> str:
        if self.mode == "vc":
            mid = (f"Vt {self.vt_ml_kg:.1f} ml/kg = {self.vt_ml_kg * ibw:.0f} ml   "
                   f"RR {self.rr:.0f}   I:E 1:{self.ie:g}")
        else:
            tag = "PS" if self.mode == "psv" else "dP"
            mid = f"{tag} {self.dp:.0f} cmH2O   RR {self.rr:.0f}   I:E 1:{self.ie:g}"
        return f"{self.label}   FiO2 {self.fio2:.2f}   PEEP {self.peep:.0f}   " + mid


@dataclass
class Patient:
    """Body habitus, respiratory mechanics, lung anatomy, metabolism, drive."""

    name: str = "adult"
    sex: str = "m"
    height_cm: float = 172.0
    weight_kg: float = 75.0
    compliance: float = 60.0             # fully recruited static compliance, mL/cmH2O
    resistance: float = 5.0              # cmH2O/L/s
    exp_limit: float = 1.0               # >1 = flow limitation, slower exhalation
    severity: float = 0.0                # 0 healthy .. 1 heavily diseased
    affected_fraction: float | None = None  # share of units that are diseased
    heterogeneity: float = 0.30          # ventilation scatter between units
    perfusion_gradient: float = 1.0      # gravity gradient of blood flow
    vascular_loss: float = 0.0           # capillary bed obliterated (PE, emphysema)
    base_shunt: float = 0.02             # anatomical venous admixture
    pfo: float = 0.0                     # patent foramen ovale: right-to-left shunt
    #                                        (opens wider as mean airway pressure loads the RV)
    lung_overrides: dict = field(default_factory=dict)  # fine-tuning of unit anatomy
    hpv: float = 0.35                    # hypoxic pulmonary vasoconstriction, 0..0.8
    n_units: int = 24
    hgb: float = 15.0                    # g/dL
    co_l_min: float = 5.0                # baseline cardiac output, L/min
    rv_sensitivity: float = 1.0          # how strongly overdistension/PEEP cuts CO
    vco2: float | None = None            # mL/min CO2 production; None -> from IBW
    rq: float = 0.8
    temp_c: float = 37.0
    hco3: float = 24.0                   # metabolic set point, before any adaptation
    renal: float = 1.0                   # 1 = normal compensation, 0 = none
    acid_mmol_h: float = 0.0             # net acid load above renal excretion
    bicarb_mmol_h: float = 0.0           # exogenous bicarbonate
    sedation: float = 1.0                # 1 = no respiratory drive, 0 = fully awake
    set_point: float = 34.0              # PaCO2 below which the CO2 drive is off
    muscle: float = 1.0                  # inspiratory muscle strength, 0 = paralysed
    ve_rest: float | None = None         # resting minute ventilation; None -> from IBW
    rr_rest: float = 15.0                # rate of the awake breathing pattern
    frc_l: float | None = None           # None -> from IBW
    notes: str = ""

    def __post_init__(self) -> None:
        self.ibw = ibw_kg(self.height_cm, self.sex)
        if self.vco2 is None:
            self.vco2 = 3.3 * self.ibw           # ~220 mL/min at 68 kg IBW
        if self.frc_l is None:
            self.frc_l = FRC_PER_KG * self.ibw
        if self.ve_rest is None:
            self.ve_rest = 0.105 * self.ibw      # ~7.1 L/min at 68 kg
        self.lung = Lung(
            n_units=self.n_units, severity=self.severity,
            affected_fraction=self.affected_fraction,
            heterogeneity=self.heterogeneity,
            perfusion_gradient=self.perfusion_gradient,
            vascular_loss=self.vascular_loss,
            base_shunt=self.base_shunt,
            **self.lung_overrides,
        )

    @property
    def weight(self) -> float:
        return self.weight_kg if self.weight_kg > 0 else self.ibw

    @property
    def buffer_space_l(self) -> float:
        return BUFFER_SPACE_FRACTION * max(self.weight_kg, 1.0)

    def describe(self) -> str:
        return (f"{self.name}: {self.sex} {self.height_cm:.0f} cm, {self.weight_kg:.0f} kg "
                f"(IBW {self.ibw:.0f} kg)\n"
                f"    C {self.compliance:.0f} ml/cmH2O   Raw {self.resistance:.0f} cmH2O/L/s"
                f"   Hb {self.hgb:.1f} g/dL   CO {self.co_l_min:.1f} L/min"
                f"   V'CO2 {self.vco2:.0f} ml/min   FRC {self.frc_l:.1f} L\n"
                f"    diseased units {100*self.lung.affected_frac:.0f}%"
                f"   HPV {self.hpv:.2f}   flow limit {self.exp_limit:.1f}"
                f"   sedation {self.sedation:.2f}"
                + (f"\n    {self.notes}" if self.notes else ""))


# ---------------------------------------------------------------------------
# result of one evaluation
# ---------------------------------------------------------------------------
@dataclass
class Result:
    t_min: float = 0.0

    # acid-base / gases
    ph: float = 7.40
    paco2: float = 40.0
    pao2: float = 90.0
    sao2: float = 0.97
    hco3: float = 24.0
    be: float = 0.0
    pvco2: float = 45.0
    pvo2: float = 40.0

    # mechanics
    mode: str = "vc"
    vt_ml: float = 500.0
    vt_ml_kg: float = 8.0
    rr: float = 12.0
    ti_s: float = 0.6
    ve_l_min: float = 6.0
    va_l_min: float = 4.5
    vd_vt: float = 0.30
    peco2: float = 30.0
    peep: float = 5.0
    auto_peep: float = 0.0
    peep_total: float = 5.0
    pplat: float = 15.0
    pip: float = 18.0
    dp_stat: float = 10.0
    pmean: float = 9.0
    peak_flow: float = 40.0
    c_dyn: float = 50.0
    local_stress: float = 12.0
    mp_j_min: float = 0.6
    mp_elas: float = 0.5
    mp_res: float = 0.1

    # oxygenation and indices
    fio2: float = 0.5
    pf: float = 100.0
    oi: float = 10.0
    shunt_pct: float = 1.0
    vq_low_pct: float = 0.0
    vq_high_pct: float = 0.0
    recruit_pct: float = 100.0
    aa_gradient: float = 10.0
    cao2: float = 0.19
    cvao2_diff: float = 0.05
    do2: float = 950.0
    co_l_min: float = 5.0
    drive: float = 0.0
    effort_ratio: float = 1.0     # delivered / wanted tidal volume in PSV

    flags: list = field(default_factory=list)
    units: list = field(default_factory=list)


# ---------------------------------------------------------------------------
# respiratory drive
# ---------------------------------------------------------------------------
def minute_target(paco2: float, sao2: float, patient: Patient) -> float:
    """Minute ventilation this patient's own drive would produce (L/min).

    A classical CO2-response line through the apnoea threshold plus a modest
    hypoxic component, both damped by `sedation`.
    """
    if patient.sedation >= 1.0:
        return 0.0
    ve = DRIVE_SLOPE * max(0.0, paco2 - patient.set_point)
    hyp = HYPOXIC_GAIN * max(0.0, 0.90 - sao2) / 0.08
    return max(0.0, (1.0 - patient.sedation) * ve * (1.0 + hyp))


def drive_index(paco2: float, sao2: float, patient: Patient) -> float:
    """Respiratory drive, normalised: 1 = resting awake ventilation, 0 = apnoea."""
    return minute_target(paco2, sao2, patient) / max(patient.ve_rest, 0.5)


# ---------------------------------------------------------------------------
# the operating point
# ---------------------------------------------------------------------------
def _evaluate(
    patient: Patient,
    vent: Vent,
    *,
    hco3: float,
    ph_for_o2: float,
    drive: float,
    t_min: float = 0.0,
    keep_units: bool = False,
    state: tuple[float, float, float] | None = None,
) -> Result:
    vent = vent.clamped()
    lung = patient.lung
    n = lung.n
    ibw = patient.ibw
    vco2 = patient.vco2

    stiff_sum = sum(lung.stiff) or 1.0
    c0 = [patient.compliance / 1000.0 * lung.stiff[i] / stiff_sum for i in range(n)]

    # ---- breath timing; spontaneous modes take their pattern from drive ---
    ve_target = drive * patient.ve_rest
    if vent.mode == "psv" and ve_target > 1e-6:
        rr_pat = _clamp(patient.rr_rest * (ve_target / patient.ve_rest) ** 0.6, 4.0, 45.0)
        rr_eff = max(rr_pat, vent.rr)          # machine backup rate
        vt_cap = ve_target / rr_pat if rr_pat > 1e-6 else None
        pmus = _clamp(PMUS_MAX * (drive - 1.0), 0.0, PMUS_CAP) * patient.muscle
        dp_eff = vent.dp + pmus
    else:
        rr_eff = vent.rr
        vt_cap = None
        dp_eff = vent.dp
    tt = 60.0 / rr_eff

    opened = [1.0] * n
    od = [0.0] * n
    stress = [0.0] * n
    shares = list(c0)
    sum_shares = sum(c0) or 1e-12
    vt_l = vent.vt_ml_kg * ibw / 1000.0
    ti = min(max(vent.ti_s or tt / (1.0 + vent.ie), 0.25), max(tt - 0.15, 0.25))
    te = tt - ti
    c_tot = patient.compliance / 1000.0
    auto = 0.0
    peep_t = vent.peep
    plat = elastic = dp_recruit = 0.0

    for _ in range(24):
        c_tot = max(sum(c0[i] * opened[i] * (1.0 - OD_STIFFEN * od[i]) for i in range(n)), 1e-5)
        tau = patient.resistance * c_tot                 # inspiratory time constant, s
        tau_exp = tau * max(patient.exp_limit, 1e-3)
        if vent.mode == "psv":
            ti = _clamp(0.45 + 1.4 * vt_l, 0.45, max(0.85 * tt, 0.3))
        else:
            ti = _clamp(vent.ti_s or tt / (1.0 + vent.ie), 0.25, max(tt - 0.15, 0.25))
        te = max(tt - ti, 1e-4)
        # dynamic hyperinflation, single time constant at steady state
        if vt_l > 0:
            auto = (vt_l / c_tot) / max(math.exp(min(te / tau_exp, 40.0)) - 1.0, 1e-9)
        else:
            auto = 0.0
        auto = min(auto, 40.0)
        peep_t = vent.peep + auto

        if vent.mode == "vc":
            vt_l = vent.vt_ml_kg * ibw / 1000.0
        else:
            vt_l = dp_eff * c_tot * (1.0 - math.exp(-ti / tau) if tau > 1e-9 else 1.0)
            if vt_cap is not None:
                vt_l = min(vt_l, vt_cap)

        elastic = vt_l / c_tot
        plat = peep_t + elastic
        dp_recruit = peep_t + RECRUIT_TIDAL * elastic

        # Recruitment is damped: in volume control a closed unit raises the very
        # pressure that would open it, so an undamped iteration can run away.
        for i in range(n):
            tgt_open = _clamp((dp_recruit - lung.p_open[i]) / lung.recruit_width, 0.0, 1.0)
            opened[i] += 0.4 * (tgt_open - opened[i])
        shares = [c0[i] * opened[i] * (1.0 - OD_STIFFEN * od[i]) * lung.vent_het[i]
                  for i in range(n)]
        sum_shares = max(sum(shares), 1e-12)
        for i in range(n):
            stress[i] = vt_l * shares[i] / sum_shares / max(c0[i], 1e-12)
            dist = peep_t * (1.0 - NONDEP_TRANSMISSION * lung.z[i])
            target = _clamp((stress[i] + dist - lung.od_limit[i]) / OD_WIDTH, 0.0, 1.0)
            if opened[i] <= 1e-9:
                target = 0.0
            od[i] += 0.5 * (target - od[i])

    # ---- perfusion: gravity, minus overdistension, minus hypoxic vasoconstriction
    mean_od = sum(lung.perf[i] * od[i] for i in range(n))
    perf = [lung.perf[i] * (1.0 - OD_CAPILLARY * od[i])
            * (1.0 - patient.hpv * (1.0 - opened[i])) for i in range(n)]
    if sum(perf) <= 1e-9:
        perf = list(lung.perf)

    # ---- ventilation distribution ----------------------------------------
    vt_dl = ANATOMIC_VD_PER_KG * ibw / 1000.0                 # anatomic dead space, L
    ve_l_min = rr_eff * vt_l
    va_l_min = max(rr_eff * max(vt_l - vt_dl, 0.0), 1e-4)
    va_ml = [va_l_min * 1000.0 * shares[i] / sum_shares for i in range(n)]

    pmean = _mean_airway_pressure(peep_t, elastic, ti, te, tt, vt_l, tau_exp,
                                  patient.resistance)
    peak_flow = (vt_l / ti * 60.0) if ti > 1e-9 else 0.0      # L/min
    resistive = patient.resistance * peak_flow / 60.0         # cmH2O
    pip = plat + resistive
    co_l_min = patient.co_l_min * (
        1.0 - patient.rv_sensitivity * (OD_CO * mean_od
                                        + OD_MEAN_CO * max(pmean - OD_MEAN_BASE, 0.0) / 10.0))
    co_l_min = max(co_l_min, 0.40 * patient.co_l_min)

    # A PFO shunts more the harder the right ventricle has to push.
    extra_shunt = 0.0
    if patient.pfo > 0.0:
        extra_shunt = patient.pfo * _clamp(1.0 + (pmean - 12.0) / 24.0, 0.0, 2.5)

    ex = exchange(
        lung, va_ml, perf,
        vco2_ml_min=vco2, co_ml_min=co_l_min * 1000.0, fio2=vent.fio2,
        hgb=patient.hgb, ph=ph_for_o2, rq=patient.rq, temp_c=patient.temp_c,
        recruited=opened, overdist=od, extra_shunt=extra_shunt,
    )

    paco2, pao2, sao2 = ex.paco2, ex.pao2, ex.sao2
    if state is not None:
        paco2, pao2, sao2 = state
    ph = phys.ph_from(hco3, paco2)
    be = phys.base_excess(hco3, ph, patient.hgb)

    # ---- indices ----------------------------------------------------------
    c_dyn = vt_l / elastic * 1000.0 if elastic > 1e-9 else 0.0
    dp_stat = elastic
    if ve_l_min > 1e-6:
        faco2 = vco2 / (ve_l_min * 1000.0 * phys.BTPS_TO_STPD)
        peco2 = faco2 * phys.PB
    else:
        peco2 = 0.0
    vd_vt = _clamp(1.0 - peco2 / paco2 if paco2 > 1e-6 else 1.0, 0.0, 0.99)
    pf = pao2 / vent.fio2 if vent.fio2 > 0.05 else 0.0
    oi = vent.fio2 * pmean * 100.0 / pao2 if pao2 > 1.0 else 9999.0
    recruit_pct = 100.0 * sum(lung.perf[i] * opened[i] for i in range(n))
    max_stress = max((stress[i] * opened[i] + peep_t * (1.0 - NONDEP_TRANSMISSION * lung.z[i])
                      for i in range(n)), default=0.0)
    mp_elas = rr_eff * 0.098 * 0.5 * elastic * vt_l
    mp_res = rr_eff * 0.098 * patient.resistance * vt_l * vt_l / max(ti, 1e-6)
    mp_j_min = mp_elas + mp_res
    cao2 = phys.o2_content(pao2, patient.hgb, ph, paco2, patient.temp_c)
    cvo2 = max(cao2 - (vco2 / patient.rq) / max(co_l_min * 1000.0, 1.0), 0.02)
    do2 = co_l_min * 1000.0 * cao2
    pao2_ideal = phys.alveolar_gas(phys.pio2(vent.fio2), paco2, patient.rq, vent.fio2)
    aa = max(pao2_ideal - pao2, 0.0)

    return Result(
        t_min=t_min, ph=ph, paco2=paco2, pao2=pao2, sao2=sao2, hco3=hco3, be=be,
        pvco2=ex.pvco2, pvo2=ex.pvo2,
        mode=vent.mode, vt_ml=vt_l * 1000.0, vt_ml_kg=vt_l * 1000.0 / ibw, rr=rr_eff,
        ti_s=ti, ve_l_min=ve_l_min, va_l_min=va_l_min, vd_vt=vd_vt, peco2=peco2,
        peep=vent.peep, auto_peep=auto, peep_total=peep_t, pplat=plat, pip=pip,
        dp_stat=dp_stat, pmean=pmean, peak_flow=peak_flow, c_dyn=c_dyn,
        local_stress=max_stress, mp_j_min=mp_j_min, mp_elas=mp_elas, mp_res=mp_res,
        fio2=vent.fio2, pf=pf, oi=oi, shunt_pct=100.0 * ex.shunt_frac,
        vq_low_pct=100.0 * ex.vq_low, vq_high_pct=100.0 * ex.vq_high,
        recruit_pct=recruit_pct, aa_gradient=aa, cao2=cao2,
        cvao2_diff=cao2 - cvo2, do2=do2, co_l_min=co_l_min, drive=drive,
        effort_ratio=(vt_l / vt_cap if vt_cap and vt_cap > 1e-6 else 1.0),
        units=ex.units if keep_units else [],
    )


def solve_steady(
    patient: Patient,
    vent: Vent,
    *,
    hco3: float | None = None,
    paco2_state: float = 40.0,
    sao2_state: float = 0.97,
    drive: float | None = None,
    t_min: float = 0.0,
    keep_units: bool = False,
    state: tuple[float, float, float] | None = None,
) -> Result:
    """Evaluate one operating point of patient + ventilator.

    `hco3` is the current metabolic bicarbonate, which with `paco2_state` fixes
    the pH used for the O2 dissociation curve (a small Bohr effect, but real).

    With `state=(PaCO2, PaO2, SaO2)` the mechanics and gas exchange are solved
    as usual but *those* values are reported, so a time-stepped simulation can
    show a real blood gas while it is still on its way to equilibrium.

    In pressure support the patient's drive and the resulting blood gas are
    coupled: drive sets the ventilation, the ventilation sets the gas, the gas
    sets the drive.  Pass `drive` to pin it (the simulation does, per step) or
    leave it `None` to solve the self-consistent value.
    """
    vent = vent.clamped()
    hco3 = patient.hco3 if hco3 is None else hco3
    ph_for_o2 = phys.ph_from(hco3, max(paco2_state, 5.0))

    if vent.mode != "psv" or drive is not None:
        d = drive_index(paco2_state, sao2_state, patient) if drive is None else drive
        res = _evaluate(patient, vent, hco3=hco3, ph_for_o2=ph_for_o2, drive=d,
                        t_min=t_min, keep_units=keep_units, state=state)
        res.flags = result_flags(res, patient)
        return res
    res = _drive_equilibrium(patient, vent, hco3, ph_for_o2, drive_index(
        paco2_state, sao2_state, patient), t_min=t_min, keep_units=keep_units, state=state)
    res.flags = result_flags(res, patient)
    return res


def _drive_equilibrium(patient: Patient, vent: Vent, hco3: float, ph_for_o2: float,
                       seed: float, *, t_min: float = 0.0, keep_units: bool = False,
                       state: tuple[float, float, float] | None = None) -> Result:
    """Solve drive <-> PaCO2 self-consistently.

    A naive fixed-point iteration does not work here: 'drive -> resulting
    PaCO2' is steeply decreasing and the CO2 response line is steeply
    increasing, so the loop gain is far above 1 and the iteration oscillates
    instead of settling.  The *difference* between the two, however, is
    monotone increasing in drive, so it has exactly one root and bisecting for
    it is robust.  Cost is ~30 evaluations of the operating point.
    """
    def gap(d: float) -> tuple[float, Result]:
        r = _evaluate(patient, vent, hco3=hco3, ph_for_o2=ph_for_o2, drive=d,
                      t_min=t_min, keep_units=keep_units, state=state)
        return d - drive_index(r.paco2, r.sao2, patient), r

    g0, res0 = gap(0.0)
    if g0 >= 0.0:                 # the machine alone keeps the CO2 under the
        res0.drive = 0.0          # apnoea threshold: no drive is needed at all
        return res0

    lo, res_lo, hi, res_hi = 0.0, res0, None, None
    d = max(seed, 0.5)
    for _ in range(8):
        gd, rd = gap(d)
        if gd >= 0.0:
            hi, res_hi = d, rd
            break
        lo, res_lo = d, rd
        d *= 2.0
    if hi is None:                # never crosses: the patient's own limit is
        res_hi = rd               # not enough to satisfy the drive they have
        res_hi.drive = d
        return res_hi

    for _ in range(20):
        mid = 0.5 * (lo + hi)
        gm, rm = gap(mid)
        if gm >= 0.0:
            hi, res_hi = mid, rm
        else:
            lo, res_lo = mid, rm
    res_hi.drive = hi
    return res_hi


def _mean_airway_pressure(peep_t, elastic, ti, te, tt, vt_l, tau_exp, resistance) -> float:
    """Analytic mean airway pressure over one breath (constant inspiratory flow)."""
    if tt <= 1e-9:
        return peep_t
    flow = vt_l / ti if ti > 1e-9 else 0.0
    resistive = resistance * flow
    mean_insp = peep_t + resistive + 0.5 * elastic
    if te > 1e-9 and tau_exp > 1e-9:
        frac = (tau_exp / te) * (1.0 - math.exp(-min(te / tau_exp, 40.0)))
    else:
        frac = 1.0
    mean_exp = peep_t + elastic * frac
    return (ti * mean_insp + te * mean_exp) / tt


# ---------------------------------------------------------------------------
# clinical flags
# ---------------------------------------------------------------------------
def result_flags(res: Result, patient: Patient) -> list[tuple[str, str]]:
    f: list[tuple[str, str]] = []
    add = f.append
    if res.pplat > 30:
        add(("alert", f"Pplat {res.pplat:.0f} > 30 cmH2O"))
    if res.dp_stat > 17:
        add(("alert", f"driving pressure {res.dp_stat:.0f} > 17 cmH2O"))
    elif res.dp_stat > 15:
        add(("warn", f"driving pressure {res.dp_stat:.0f} > 15 cmH2O"))
    if res.vt_ml_kg > 9.5:
        add(("warn", f"Vt {res.vt_ml_kg:.1f} ml/kg IBW is large"))
    if res.local_stress > 28:
        add(("warn", f"peak local stress {res.local_stress:.0f} cmH2O (regional overstretch)"))
    if res.auto_peep > 8:
        add(("alert", f"auto-PEEP {res.auto_peep:.0f} cmH2O: breath stacking"))
    elif res.auto_peep > 4:
        add(("warn", f"auto-PEEP {res.auto_peep:.0f} cmH2O"))
    if res.mp_j_min > 10:
        add(("warn", f"tidal mechanical power {res.mp_j_min:.1f} J/min"))
    if res.sao2 < 0.90:
        add(("alert", f"SaO2 {100*res.sao2:.0f}% - hypoxaemia"))
    if res.pao2 > 130 and res.fio2 > 0.35:
        add(("warn", f"PaO2 {res.pao2:.0f} on FiO2 {res.fio2:.2f}: hyperoxia, wean FiO2"))
    if res.ph < 7.30:
        add(("alert", f"pH {res.ph:.2f} acidemia"))
    elif res.ph > 7.50:
        add(("alert", f"pH {res.ph:.2f} alkalemia"))
    if res.paco2 > 55 and res.hco3 < 30:
        add(("warn", f"CO2 {res.paco2:.0f} with little bicarbonate buffer"))
    if res.vd_vt > 0.65:
        add(("warn", f"Vd/Vt {res.vd_vt:.2f}: high dead space"))
    if res.co_l_min < 0.75 * patient.co_l_min:
        add(("warn", f"cardiac output {res.co_l_min:.1f} L/min, down "
                     f"{100*(1-res.co_l_min/patient.co_l_min):.0f}% from pressure effects"))
    if res.drive > 1.8:
        add(("warn", f"respiratory drive {res.drive:.1f}x resting: high work of breathing"))
    if res.mode == "psv" and res.effort_ratio < 0.97:
        add(("warn", f"support-limited: only {100*res.effort_ratio:.0f}% of the "
                     f"targeted tidal volume is delivered"))
    if patient.severity >= 0.25 and res.vt_ml_kg > 8.0:
        add(("warn", f"stiff lung at {res.vt_ml_kg:.1f} ml/kg - protective ceiling is 8"))
    if (res.sao2 >= 0.90 and res.pao2 >= 62 and 7.35 <= res.ph <= 7.48
            and res.dp_stat <= 15 and res.auto_peep <= 4):
        add(("ok", "oxygenation, ventilation and stress targets all met"))
    return f
