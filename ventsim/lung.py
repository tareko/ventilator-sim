"""Multi-compartment V/Q lung.

Structure
---------
The lung is a stack of parallel units ordered non-dependent (i = 0) to
dependent (i = N-1), as West's zones are drawn in every textbook.  A `severity`
parameter turns a dependent *fraction* of the units into diseased ones:

  healthy units    closing pressure ~1-3 cmH2O (open at any sane PEEP),
                   normal compliance, overdistend at a *low* stress
  diseased units   closing pressure 7-30 cmH2O (recruit gradually with PEEP),
                   stiff, and protected from overdistension by the surrounding
                   oedema/consolidation - which is exactly the ARDS picture of
                   vulnerable non-dependent lung beside protected dependent
                   collapse.

Ventilation is distributed by local compliance and perfusion by the gravity
gradient, so low-V/Q dependent units and high-V/Q overdistended units fall out
of the geometry instead of being bolted on as a fake "A-a gradient".

Gas exchange
------------
CO2 (linear blood curve -> closed form isopleth):

      Pc_i = PvCO2 / (1 + Vstpd_i / (Q_i * beta * (Pb - PH2O)))

O2 - a real mass balance per unit, not the alveolar gas equation:

      Q_i * (CcO2_i - CvO2) = Vstpd_i * (FIO2 - PAO2_i / (Pb - PH2O))

The difference matters.  With the gas equation, a unit whose V/Q approaches
zero is handed a *negative* PO2 and the blood leaving it carries no oxygen at
all; with the mass balance the same unit simply cannot deliver much oxygen, so
its blood ends up venous - which is what a poorly ventilated unit really does.
Solved by bisection (the left-hand side is monotone in PO2), with mixed venous
content iterated because shunt blood feeds back into it.

Because the whole thing is derived from mass balance, anaemia, low cardiac
output and hypoxic vasoconstriction all change the blood gas the way they do at
the bedside, with no extra fudge factors.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from . import physiology as phys

MAX_PACO2 = 160.0     # numerical ceiling when alveolar ventilation approaches zero


def _clamp(x: float, lo: float, hi: float) -> float:
    return lo if x < lo else hi if x > hi else x


@dataclass
class Lung:
    """Anatomy of one patient's lung.  Unit tables are built once."""

    n_units: int = 24
    severity: float = 0.0              # 0 healthy .. 1 heavily diseased
    affected_fraction: float | None = None  # defaults to 0.55*severity
    heterogeneity: float = 0.30        # ventilation scatter between units
    perfusion_gradient: float = 1.0    # gravity gradient of blood flow
    vascular_loss: float = 0.0         # share of the capillary bed obliterated
    #                                        (pulmonary embolism, emphysema)
    healthy_closing: float = 3.0       # closing pressure of normal units, cmH2O
    affected_closing_lo: float = 7.0   # diseased units: recruit between ...
    affected_closing_hi: float = 30.0  # ... and here
    od_lo_healthy: float = 17.0        # stress limit of normal units
    od_span_healthy: float = 14.0
    od_lo_affected: float = 30.0       # diseased units tolerate much more
    od_span_affected: float = 25.0
    stiff_affected: float = 0.7        # compliance of diseased units, relative
    recruit_width: float = 5.0         # pressure range over which a unit opens
    base_shunt: float = 0.02           # anatomical + bronchial venous admixture
    max_closing: float | None = None   # optional hard override of closing spread

    def __post_init__(self) -> None:
        self.n_units = max(4, int(self.n_units))
        sev = _clamp(self.severity, 0.0, 1.0)
        self.affected_frac = (0.55 * sev) if self.affected_fraction is None \
            else _clamp(self.affected_fraction, 0.0, 1.0)
        self.build()

    def build(self) -> None:
        n = self.n_units
        af = self.affected_frac
        edge = max(1.5 / n, 0.02)          # soft boundary between healthy and diseased
        start = 1.0 - af
        self.z = [0.0] * n
        self.perf = [0.0] * n
        self.p_open = [0.0] * n
        self.od_limit = [0.0] * n
        self.stiff = [1.0] * n
        self.vent_het = [0.0] * n
        self.aff = [0.0] * n

        raw = [0.0] * n
        self.perf_scale = [1.0] * n
        tot = 0.0
        for i in range(n):
            z = (i + 0.5) / n
            self.z[i] = z
            # capillary loss lands on a deterministic scatter of units, so some
            # are well ventilated and badly perfused (dead space) by design
            mask = 0.5 + 0.5 * math.sin((i + 1) * 1.7138750)
            scale = _clamp(1.0 - 2.0 * self.vascular_loss * mask, 0.02, 1.0)
            self.perf_scale[i] = scale
            p = (0.55 + self.perfusion_gradient * z) * scale
            raw[i] = p
            tot += p
            aff = _clamp((z - start) / edge, 0.0, 1.0)
            self.aff[i] = aff

            healthy_close = self.healthy_closing * (0.25 + 0.75 * z)
            span = max(af, 1e-6)
            dis_close = (self.affected_closing_lo
                         + (self.affected_closing_hi - self.affected_closing_lo)
                         * _clamp((z - start) / span, 0.0, 1.0))
            if self.max_closing is not None:
                dis_close = self.max_closing * z
            self.p_open[i] = healthy_close + aff * (dis_close - healthy_close)

            od_h = self.od_lo_healthy + self.od_span_healthy * z
            od_a = self.od_lo_affected + self.od_span_affected * z
            self.od_limit[i] = od_h + aff * (od_a - od_h)
            self.stiff[i] = 1.0 + aff * (self.stiff_affected - 1.0)
            # golden-angle scatter: deterministic, identical value in JS
            self.vent_het[i] = 1.0 + self.heterogeneity * math.sin((i + 1) * 2.3999632)

        for i in range(n):
            self.perf[i] = raw[i] / tot

    @property
    def n(self) -> int:
        return self.n_units


@dataclass
class UnitRow:
    z: float
    perf: float
    va: float
    vq: float
    pc_co2: float
    po2: float
    cc_o2: float
    shunt: bool
    diseased: float = 0.0
    recruited: float = 1.0
    overdist: float = 0.0


@dataclass
class ExchangeResult:
    pao2: float
    paco2: float
    sao2: float
    cao2: float
    cvo2: float
    pvco2: float
    pvo2: float
    vo2_ml_min: float
    shunt_frac: float     # perfusion share of unventilated units (+ anatomical)
    vq_high: float        # ventilation share above V/Q 2  (dead-space-like)
    vq_low: float         # ventilation share below V/Q 0.6
    units: list = field(default_factory=list)


def _solve_unit_po2(a: float, cvo2: float, fio2: float, b_cap: float,
                    p50_mm: float, dry: float, alpha: float) -> float:
    """PO2 of one unit from its O2 mass balance, by bisection.

      B*so2(p) + alpha*p = cvo2 + a*(FIO2 - p/dry)

    so2(p) is monotone increasing while the right-hand side minus the left is
    monotone decreasing, so the zero crossing is unique and bisection is safe.
    """
    k = alpha + a / dry
    c0 = cvo2 + a * fio2
    hi = fio2 * dry
    if hi <= 0.0:
        return 0.0

    def g(p: float) -> float:
        s_rhs = (c0 - k * p) / b_cap
        if p <= 0.0:
            return -s_rhs
        ln_x = phys.HILL_N * math.log(p / p50_mm)
        if ln_x > 40.0:
            s_h = 1.0
        elif ln_x < -40.0:
            s_h = math.exp(ln_x)
        else:
            e = math.exp(ln_x)
            s_h = e / (1.0 + e)
        return s_h - s_rhs

    lo_p, hi_p = 0.0, hi
    if g(hi_p) < 0.0:            # only reachable with pathological parameters
        return hi_p
    for _ in range(26):
        mid = 0.5 * (lo_p + hi_p)
        if g(mid) < 0.0:
            lo_p = mid
        else:
            hi_p = mid
    return 0.5 * (lo_p + hi_p)


def exchange(
    lung: Lung,
    va_ml: list[float],
    perf: list[float],
    *,
    vco2_ml_min: float,
    co_ml_min: float,
    fio2: float,
    hgb: float = 15.0,
    ph: float = 7.40,
    rq: float = 0.8,
    temp_c: float = 37.0,
    recruited: list[float] | None = None,
    overdist: list[float] | None = None,
    extra_shunt: float = 0.0,
) -> ExchangeResult:
    """Solve the compartment system.

    ``va_ml``: alveolar ventilation per unit, mL/min BTPS.
    ``perf``: perfusion distribution (any scale, it is normalised).
    ``extra_shunt``: intracardiac / intrapulmonary right-to-left flow on top of
    the anatomical admixture (a PFO opened up by a pressure-loaded RV).
    """
    n = lung.n
    if len(va_ml) != n or len(perf) != n:
        raise ValueError("per-unit vectors must match lung.n_units")

    tot_perf = sum(perf)
    if tot_perf <= 0:
        raise ValueError("perfusion distribution is empty")
    w = [p / tot_perf for p in perf]

    # Anatomical admixture: bronchial/Thebesian flow (plus any intracardiac
    # shunt handed in) that never sees gas.
    fs = _clamp(lung.base_shunt + extra_shunt, 0.0, 0.5)
    w_eff = [p * (1.0 - fs) for p in w]

    co = max(co_ml_min, 500.0)
    vo2_ml_min = vco2_ml_min / rq
    k_o2 = vo2_ml_min / co                      # arterio-venous O2 content difference
    beta = phys.BETA_CO2
    dry = phys.DRY
    denom = beta * dry

    vst = [max(v, 0.0) * phys.BTPS_TO_STPD for v in va_ml]

    # ---- CO2, closed form ------------------------------------------------
    # Pc_i = PvCO2/(1+a_i) and PvCO2 = PaCO2 + K give PaCO2 = S*(PaCO2 + K).
    a_list = [0.0] * n
    q_list = [0.0] * n
    shunt_i = [False] * n
    s_sum = fs
    shunt = fs
    for i in range(n):
        q = w_eff[i] * co
        q_list[i] = q
        a = (vst[i] / (q * denom)) if q > 1e-9 else float("inf")
        a_list[i] = a
        s_sum += w_eff[i] / (1.0 + a)
        # Blood from a unit whose V/Q is under ~0.05 leaves essentially
        # untouched, whether or not its airway is technically open.
        if q > 1e-9 and va_ml[i] / q < 0.05:
            shunt_i[i] = True
            shunt += w_eff[i]
    k_co2 = vco2_ml_min / (co * beta)
    if s_sum >= 0.999:
        paco2 = MAX_PACO2
    else:
        paco2 = min(s_sum * k_co2 / (1.0 - s_sum), MAX_PACO2)
    pvco2 = paco2 + k_co2

    # ---- O2: mass balance per unit, venous content iterated --------------
    b_cap = phys.HUNTER * (hgb / 100.0)
    alpha = phys.ALPHA_O2
    p50_mm = phys.p50(ph, pvco2, temp_c)
    vented = [q_list[i] > 1e-9 and va_ml[i] > 1e-9 and not shunt_i[i]
              for i in range(n)]
    cvo2 = max(min(0.13, b_cap * 0.65), 0.02)   # mixed-venous starting guess
    p_list = [0.0] * n
    cao2 = 0.15
    for _ in range(12):
        # anatomical admixture and unventilated units both carry venous blood
        cao2 = shunt * cvo2
        for i in range(n):
            if not vented[i]:
                continue
            p_list[i] = _solve_unit_po2(vst[i] / q_list[i], cvo2, fio2,
                                        b_cap, p50_mm, dry, alpha)
            cao2 += w_eff[i] * phys.o2_content(p_list[i], hgb, ph, pvco2, temp_c)
        cvo2_next = max(cao2 - k_o2, 0.01)
        if abs(cvo2_next - cvo2) < 1e-7:
            cvo2 = cvo2_next
            break
        cvo2 = cvo2 + 0.55 * (cvo2_next - cvo2)

    pvo2 = phys.content_to_o2(cvo2, hgb, ph, pco2=pvco2, temp_c=temp_c)[0]
    for i in range(n):
        if not vented[i]:
            # unventilated units leave blood venous; unperfused ones hold inspired gas
            p_list[i] = pvo2 if q_list[i] > 1e-9 else fio2 * dry

    pao2, sao2 = phys.content_to_o2(cao2, hgb, ph, pco2=paco2, temp_c=temp_c)

    # ---- per-unit report rows -------------------------------------------
    vent_total = sum(va_ml)
    vq_high = vq_low = 0.0
    rows: list[UnitRow] = []
    for i in range(n):
        q = q_list[i]
        vq = (va_ml[i] / q) if q > 1e-9 else float("inf")
        if va_ml[i] > 1e-9 and vent_total > 0:
            frac = va_ml[i] / vent_total
            if vq > 2.0:
                vq_high += frac
            elif vq < 0.6:
                vq_low += frac
        rows.append(UnitRow(
            z=lung.z[i], perf=w_eff[i], va=va_ml[i], vq=vq,
            pc_co2=0.0 if a_list[i] == float("inf") else pvco2 / (1.0 + a_list[i]),
            po2=p_list[i],
            cc_o2=phys.o2_content(p_list[i], hgb, ph, pvco2, temp_c),
            shunt=shunt_i[i],
            diseased=lung.aff[i],
            recruited=1.0 if recruited is None else recruited[i],
            overdist=0.0 if overdist is None else overdist[i],
        ))

    return ExchangeResult(
        pao2=pao2, paco2=paco2, sao2=sao2, cao2=cao2, cvo2=cvo2,
        pvco2=pvco2, pvo2=pvo2, vo2_ml_min=vo2_ml_min, shunt_frac=shunt,
        vq_high=vq_high, vq_low=vq_low, units=rows,
    )
