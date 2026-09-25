"""Time-stepped simulation: how long a change takes to show up in the blood gas.

All the relaxations are exact exponential integrations of first-order
compartments, so they are stable at any step size:

  PaCO2   lung/blood store (tau ~ FRC / alveolar ventilation)
          + tissue CO2 store (tau ~ 30 / cardiac output)
  PaO2    same structure but faster - the O2 store is much smaller
  HCO3    acute cellular buffering (tau ~ 0.5 min)
          + renal compensation (tau ~ 12 h, switched off by renal failure)
          + exogenous acid or bicarbonate load
  drive   respiratory drive follows the blood gases with a short neural lag

Consequence worth playing with: fix a chronically high CO2 and the bicarbonate
the kidneys grew over two days does not evaporate, so pH swings alkalotic -
post-hypercapnic alkalosis, for free.
"""

from __future__ import annotations

import math
from dataclasses import replace

from . import physiology as phys
from .model import (
    ACUTE_HCO3_PER_MMHG,
    CHRONIC_HCO3_PER_MMHG,
    CHRONIC_HYPO_PER_MMHG,
    TAU_ACUTE_MIN,
    TAU_RENAL_MIN,
    Patient,
    Result,
    Vent,
    _clamp,
    drive_index,
    result_flags,
    solve_steady,
)

MAX_STEPS = 400_000
TAU_DRIVE_MIN = 0.4           # neural lag of the respiratory centre
MAX_DRIVE_STEP_MIN = 0.1      # sub-step limit while the drive loop is closed


def _relax(x: float, target: float, tau: float, dt: float) -> float:
    """Exact first-order relaxation over dt (unconditionally stable)."""
    if tau <= 1e-9:
        return target
    return target + (x - target) * math.exp(-dt / tau)


class Simulation:
    """One patient on one ventilator, advanced in time."""

    def __init__(self, patient: Patient, vent: Vent | None = None, *,
                 hco3: float | None = None, t_min: float = 0.0,
                 paco2: float | None = None, pao2: float | None = None):
        self.patient = patient
        self.vent = (vent or Vent()).clamped()
        self._hco3_base = float(patient.hco3 if hco3 is None else hco3)
        self._hco3_buf = 0.0
        self.t_min = t_min
        base = solve_steady(patient, self.vent, hco3=self.hco3)
        self.paco2 = base.paco2 if paco2 is None else paco2
        self.pao2 = base.pao2 if pao2 is None else pao2
        self._fast_co2 = self.paco2
        self._slow_co2 = self.paco2
        self._fast_o2 = self.pao2
        self._slow_o2 = self.pao2
        self.drive = drive_index(self.paco2, self.sao2_current(), patient)

    # ---------------------------------------------------------------- info
    @property
    def hco3(self) -> float:
        """Total bicarbonate: the metabolic baseline plus acute buffering."""
        return self._hco3_base + self._hco3_buf

    @hco3.setter
    def hco3(self, value: float) -> None:
        self._hco3_base = float(value)
        self._hco3_buf = 0.0

    def sao2_current(self) -> float:
        ph = phys.ph_from(self.hco3, max(self.paco2, 4.0))
        return phys.so2(self.pao2, ph=ph, pco2=max(self.paco2, 4.0),
                        temp_c=self.patient.temp_c)

    def target(self, keep_units: bool = False) -> Result:
        """Where the gases are heading if nothing else changes."""
        return solve_steady(self.patient, self.vent, hco3=self.hco3,
                            paco2_state=self.paco2, sao2_state=self.sao2_current(),
                            drive=self.drive, keep_units=keep_units)

    def current(self, keep_units: bool = False) -> Result:
        """The blood gas and the mechanics right now."""
        return solve_steady(
            self.patient, self.vent, hco3=self.hco3,
            paco2_state=self.paco2, sao2_state=self.sao2_current(),
            drive=self.drive, t_min=self.t_min, keep_units=keep_units,
            state=(self.paco2, self.pao2, self.sao2_current()),
        )

    # ------------------------------------------------------------ stepping
    def step(self, dt_min: float) -> Result:
        dt_total = max(min(dt_min, 30.0), 0.0)
        if dt_total <= 0:
            return self.current()
        # Pressure support is a feedback loop: the drive reacts to the blood
        # gases and the gases react to the ventilation that drive produces.
        # Solved in steps longer than about a tenth of the drive time constant
        # the loop overshoots and the simulation starts breathing in a limit
        # cycle, so sub-step whenever that loop is actually closed.
        n = 1
        if self.vent.mode == "psv" and self.patient.sedation < 1.0:
            n = math.ceil(dt_total / MAX_DRIVE_STEP_MIN)
        h = dt_total / n
        for _ in range(n):
            self._step_once(h)
        return self.current()

    def _step_once(self, dt: float) -> Result:
        p = self.patient
        tgt = self.target()
        va = max(tgt.va_l_min, 0.15)
        co = max(tgt.co_l_min, 0.4)

        tau_fast_co2 = max(0.05, 0.6 * p.frc_l / va)
        tau_slow_co2 = 30.0 / co
        tau_fast_o2 = max(0.03, 0.5 * p.frc_l / va)
        tau_slow_o2 = 10.0 / co

        self._fast_co2 = _relax(self._fast_co2, tgt.paco2, tau_fast_co2, dt)
        self._slow_co2 = _relax(self._slow_co2, self._fast_co2, tau_slow_co2, dt)
        self.paco2 = max(0.5 * (self._fast_co2 + self._slow_co2), 3.0)

        self._fast_o2 = _relax(self._fast_o2, tgt.pao2, tau_fast_o2, dt)
        self._slow_o2 = _relax(self._slow_o2, self._fast_o2, tau_slow_o2, dt)
        self.pao2 = max(0.5 * (self._fast_o2 + self._slow_o2), 4.0)

        # ---- bicarbonate: acute cellular buffering + slow renal compensation
        self._hco3_buf = _relax(self._hco3_buf,
                                ACUTE_HCO3_PER_MMHG * (self.paco2 - 40.0), TAU_ACUTE_MIN, dt)
        if p.renal > 1e-6:
            # Compensation is asymmetric, as it is in real kidneys: retaining
            # bicarbonate against chronic hypercapnia is easy, getting rid of
            # enough to match a low PaCO2 is not.
            gain = CHRONIC_HCO3_PER_MMHG if self.paco2 >= 40.0 else CHRONIC_HYPO_PER_MMHG
            base_target = p.hco3 + gain * (self.paco2 - 40.0)
            self._hco3_base = _relax(self._hco3_base, base_target,
                                     TAU_RENAL_MIN / max(p.renal, 1e-6), dt)
        load = (p.bicarb_mmol_h - p.acid_mmol_h) / 60.0 / p.buffer_space_l  # mmol/L/min
        self._hco3_base = _clamp(self._hco3_base + load * dt, 5.0, 60.0)

        # ---- respiratory drive follows the blood gases
        drive_now = drive_index(self.paco2, self.sao2_current(), p)
        self.drive = _relax(self.drive, drive_now, TAU_DRIVE_MIN, dt)
        self.t_min += dt

    # ---------------------------------------------------------------- runs
    def run(self, minutes: float, sample_every: float = 1.0,
            dt: float | None = None) -> list[Result]:
        """Advance `minutes`; return samples every `sample_every` minutes."""
        if minutes <= 0:
            return []
        if dt is None:
            dt = min(sample_every / 4.0, 0.25)
        dt = max(dt, 1e-4)
        n_steps = max(int(round(minutes / dt)), 1)
        if n_steps > MAX_STEPS:
            raise ValueError(f"run too long: {n_steps} steps (max {MAX_STEPS}); "
                             "pass a coarser dt")
        n_samp = max(int(round(sample_every / dt)), 1)
        out: list[Result] = []
        for k in range(1, n_steps + 1):
            self.step(dt)
            if k % n_samp == 0:
                out.append(self.current())
        if not out or n_steps % n_samp:
            out.append(self.current())
        return out

    def adapt(self, hours: float) -> list[Result]:
        """Let the kidneys adapt to the current settings (chronic baseline)."""
        minutes = hours * 60.0
        return self.run(minutes, sample_every=max(minutes / 4.0, 5.0), dt=2.0)

    # --------------------------------------------------------------- editing
    def set_vent(self, **kw) -> Vent:
        self.vent = replace(self.vent, **kw).clamped()
        return self.vent

    def set_patient(self, **kw) -> Patient:
        """Replace the patient (or its parameters).  Gas state is kept."""
        self.patient = replace(self.patient, **kw)
        return self.patient

    def equilibrate(self) -> Result:
        """Jump straight to equilibrium for the current settings."""
        tgt = self.target()
        self.paco2, self.pao2 = tgt.paco2, tgt.pao2
        self._fast_co2 = self._slow_co2 = tgt.paco2
        self._fast_o2 = self._slow_o2 = tgt.pao2
        self.drive = drive_index(self.paco2, self.sao2_current(), self.patient)
        return self.current()

    def flags_now(self) -> list[tuple[str, str]]:
        return result_flags(self.current(), self.patient)
