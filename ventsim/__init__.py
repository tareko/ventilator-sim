"""vent-sim: bedside ventilator and blood-gas simulator.

Change a setting, watch the blood gases move and - just as importantly - watch
*how long* the change takes.

    from ventsim import Patient, Vent, Simulation
    p = Patient(name="ARDS moderate", severity=0.55, compliance=28)
    s = Simulation(p, Vent(mode="vc", fio2=0.9, peep=15, vt_ml_kg=6, rr=28))
    print(s.current().pao2, s.current().paco2)
"""

from .dynamics import Simulation
from .lung import Lung, UnitRow
from .model import (
    MODES,
    Patient,
    Result,
    Vent,
    drive_index,
    ibw_kg,
    minute_target,
    result_flags,
    solve_steady,
)

__all__ = [
    "Lung", "UnitRow", "Patient", "Vent", "Result", "Simulation",
    "MODES", "ibw_kg", "drive_index", "minute_target", "solve_steady", "result_flags",
]

__version__ = "0.1.0"
