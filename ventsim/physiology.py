"""Gas physiology primitives: O2/Hb dissociation, CO2 content, Henderson-Hasselbalch.

Every function is pure. Units are stated in the docstring; contents are always
mL of gas per **mL** of whole blood (so 0.20 == 20 mL/dL).

Curve choices (see README "Where the numbers come from"):

* O2-Hb dissociation: Hill equation, n = 2.7, P50 = 26.6 mmHg.  Matches the
  classic curve within ~1% over the clinical range (pO2 40 -> 75%, 60 -> 90%,
  100 -> 97%).
* Bohr / CO2 / temperature shifts applied to P50.
* CO2 content: single linear curve through the origin,
  0.0121 mL CO2 per mL blood per mmHg (48.5 mL/dL at 40 mmHg).  The CO2 curve
  is nearly linear in the physiological range and a linear curve lets the
  V/Q solver be done in closed form (Rahn isopleths).
* pH/HCO3: Henderson-Hasselbalch with pK 6.101 and alpha 0.0301 mmol/L/mmHg.
* Base excess: Van Slyke / Siggaard-Andersen, haemoglobin corrected.
"""

from __future__ import annotations

import math

# ---------------------------------------------------------------------------
# constants
# ---------------------------------------------------------------------------
PB = 760.0          # barometric pressure, mmHg
PH2O = 47.0         # water vapour pressure at 37 C, mmHg
DRY = PB - PH2O     # 713 mmHg: (PB - PH2O)

# Multiply a BTPS gas volume by this to get STPD: (PB-PH2O)/PB * 273.15/310.35
BTPS_TO_STPD = (DRY / PB) * (273.15 / 310.35)   # ~0.8257

ALPHA_O2 = 0.00003  # mL O2 dissolved per mL blood per mmHg  (== 0.003 mL/dL/mmHg)
HUNTER = 1.34       # mL O2 carried per g Hb when fully saturated
BETA_CO2 = 0.0121   # mL CO2 per mL blood per mmHg (linear whole-blood curve)

PK_CO2 = 6.101      # apparent pK of the bicarbonate buffer at 37 C
ALPHA_CO2 = 0.0301  # mmol CO2 dissolved per L plasma per mmHg

HILL_N = 2.7
P50_NORMAL = 26.6   # mmHg at pH 7.40, pCO2 40, 37 C


# ---------------------------------------------------------------------------
# pH / bicarbonate / base excess
# ---------------------------------------------------------------------------
def ph_from(hco3: float, pco2: float) -> float:
    """Henderson-Hasselbalch pH from bicarbonate (mmol/L) and PCO2 (mmHg)."""
    hco3 = max(hco3, 1e-4)
    pco2 = max(pco2, 1e-4)
    return PK_CO2 + math.log10(hco3 / (ALPHA_CO2 * pco2))


def hco3_from(ph: float, pco2: float) -> float:
    """Bicarbonate (mmol/L) matching a pH at a given PCO2."""
    return (10.0 ** (ph - PK_CO2)) * ALPHA_CO2 * pco2


def base_excess(hco3: float, ph: float, hgb: float = 15.0) -> float:
    """Base excess (mmol/L), Van Slyke / Siggaard-Andersen, Hb corrected.

    BE = (1 - 0.014*Hb) * (HCO3 - 24 + (1.43*Hb + 7.7) * (pH - 7.4))
    """
    return (1.0 - 0.014 * hgb) * (hco3 - 24.0 + (1.43 * hgb + 7.7) * (ph - 7.4))


def anion_gap(na: float, cl: float, hco3: float, albumin: float = 4.0) -> float:
    """Anion gap, corrected to a normal albumin (a low albumin hides a gap)."""
    return (na - cl - hco3) + 2.5 * (4.0 - albumin)


# ---------------------------------------------------------------------------
# O2-Hb dissociation
# ---------------------------------------------------------------------------
def p50(ph: float = 7.4, pco2: float = 40.0, temp_c: float = 37.0) -> float:
    """Position of the O2-Hb curve (mmHg). Right shift = higher value.

    The Bohr effect (pH) dominates: -0.44 per unit, the in vitro whole-blood
    figure.  CO2 has a small effect of its own through carbamino compounds,
    worth about 0.6 mmHg of P50 per 10 mmHg of PCO2 once the pH shift has
    already been counted.  Temperature moves the curve about 5.7 % per degree.
    """
    return P50_NORMAL * (
        10.0 ** (-0.44 * (ph - 7.4) + 0.0010 * (pco2 - 40.0) + 0.024 * (temp_c - 37.0))
    )


def so2(po2: float, p50_mm: float | None = None, ph: float = 7.4,
        pco2: float = 40.0, temp_c: float = 37.0) -> float:
    """Hb oxygen saturation (fraction 0..1) from dissolved O2 tension.

    Without an explicit P50 the curve is shifted for the pH, CO2 and
    temperature supplied (see `p50`), which is what lets the Bohr effect
    reach the content bookkeeping.
    """
    po2 = max(po2, 0.0)
    if po2 <= 0.0:
        return 0.0
    p50_ = p50(ph, pco2, temp_c) if p50_mm is None else p50_mm
    x = po2 ** HILL_N
    return x / (p50_ ** HILL_N + x)


def o2_content(po2: float, hgb: float = 15.0, ph: float = 7.4,
               pco2: float = 40.0, temp_c: float = 37.0) -> float:
    """O2 content (mL O2 per mL blood) = Hb bound + dissolved."""
    s = so2(po2, ph=ph, pco2=pco2, temp_c=temp_c)
    return HUNTER * (hgb / 100.0) * s + ALPHA_O2 * po2


def content_to_o2(content: float, hgb: float = 15.0, ph: float = 7.4,
                  pco2: float = 40.0, temp_c: float = 37.0,
                  po2_hi: float = 3000.0) -> tuple[float, float]:
    """Invert the O2-Hb curve: (po2, so2) from a content value.

    Bisection on the (strictly increasing) content function.  Needed because a
    naive Hill inversion loses all resolution once Hb is saturated, which is
    exactly where the dissolved-O2 part of the curve lives (high FiO2, shunt
    calculations).
    """
    if content <= 0.0:
        return 0.0, 0.0

    def c(po2: float) -> float:
        return o2_content(po2, hgb, ph, pco2, temp_c)

    lo, hi = 0.0, po2_hi
    if c(hi) < content:            # supersaturated (e.g. hyperbaric FiO2 1.0)
        bound_cap = HUNTER * (hgb / 100.0)
        return max((content - bound_cap) / ALPHA_O2, hi), 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if c(mid) < content:
            lo = mid
        else:
            hi = mid
    po2 = 0.5 * (lo + hi)
    return po2, so2(po2, ph=ph, pco2=pco2, temp_c=temp_c)


def oxygen_delivery(co_l_min: float, cao2: float) -> float:
    """Oxygen delivery (mL/min) from cardiac output (L/min) and CaO2 (mL/mL)."""
    return co_l_min * 1000.0 * cao2


# ---------------------------------------------------------------------------
# CO2 content
# ---------------------------------------------------------------------------
def co2_content(pco2: float) -> float:
    """CO2 content (mL CO2 per mL blood), linear whole-blood curve."""
    return BETA_CO2 * max(pco2, 0.0)


def pco2_from_content(content: float) -> float:
    return max(content, 0.0) / BETA_CO2


# ---------------------------------------------------------------------------
# alveolar gas equation (used per V/Q unit)
# ---------------------------------------------------------------------------
def pio2(fio2: float, pb: float = PB) -> float:
    """Inspired O2 partial pressure at the alveoli."""
    return fio2 * (pb - PH2O)


def alveolar_gas(pio2_val: float, pco2: float, rq: float = 0.8,
                 fio2: float = 0.21) -> float:
    """Full alveolar gas equation: PAO2 = PIO2 - PaCO2*(1 - FIO2*(1-RQ))/(RQ)."""
    coef = (1.0 - fio2 * (1.0 - rq)) / rq
    return max(pio2_val - pco2 * coef, 0.0)


def minute_to_alveolar_pco2(va_l_min: float, vco2_ml_min: float) -> float:
    """The textbook PaCO2 = 0.863 * VCO2 / VA (equivalent to the mass balance)."""
    if va_l_min <= 1e-6:
        return 200.0
    return 0.863 * vco2_ml_min / va_l_min
