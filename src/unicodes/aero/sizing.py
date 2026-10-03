"""Preliminary weight sizing with Roskam's mission fuel-fraction method (AE 521).

Units follow Roskam: range in nmi, speed in kts, endurance in hours,
TSFC ``c_j`` in lbf/(lbf hr); weights in any consistent unit.

    mission = [
        ("engine start", 0.99), ("taxi", 0.99), ("take-off", 0.995), ("climb", 0.98),
        ("cruise", breguet_range_fraction(3000, 430, 0.34, 26)),
        ("loiter", breguet_endurance_fraction(0.75, 0.34, 26)),
        ("descent", 0.995), ("landing", 1.0), ("shutdown", 0.995),
    ]
    est = estimate_weights(W_TO=130000, mission=mission, payload=0, crew=700)
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

# Typical phase weight fractions (Roskam, Airplane Design Part I, Table 2.1, jets)
ROSKAM_JET_FRACTIONS = {
    "engine start": 0.99,
    "taxi": 0.99,
    "take-off": 0.995,
    "climb": 0.98,
    "descent": 0.99,
    "landing": 0.992,
}


def breguet_range_fraction(range_nmi, V_kts, c_j, L_D):
    """End/start weight ratio for a jet cruise segment (Breguet range)."""
    return np.exp(-range_nmi * c_j / (V_kts * L_D))


def breguet_endurance_fraction(hours, c_j, L_D):
    """End/start weight ratio for a jet loiter segment (Breguet endurance)."""
    return np.exp(-hours * c_j / L_D)


def mach_to_kts(M, speed_of_sound_fts):
    """True airspeed in knots from Mach number and speed of sound in ft/s."""
    return M * speed_of_sound_fts * 3600 / 6076


@dataclass(frozen=True)
class WeightEstimate:
    W_TO: float
    fuel_used: float
    trapped_fuel_oil: float
    fuel: float  # used + trapped
    operating_empty_tentative: float
    empty_tentative: float
    mission_fuel_fraction: float  # landing/take-off weight ratio due to fuel only


def estimate_weights(
    W_TO: float,
    mission: Sequence[tuple[str, float]],
    payload: float = 0.0,
    crew: float = 0.0,
    payload_drop_after: str | None = None,
    trapped_fraction: float = 0.01,
) -> WeightEstimate:
    """Tentative empty weight for a guessed take-off weight.

    ``mission`` is an ordered list of ``(phase, W_end/W_start)``. If
    ``payload_drop_after`` names a phase, the payload is released at the end
    of it (as in the Hermes launch mission) and later fuel is burned at the
    lighter weight. Compare ``empty_tentative`` with a statistical empty
    weight and iterate ``W_TO`` until they match.
    """
    W = W_TO
    fuel_used = 0.0
    payload_on_board = payload
    for phase, fraction in mission:
        burned = W * (1 - fraction)
        fuel_used += burned
        W -= burned
        if phase == payload_drop_after:
            W -= payload_on_board
            payload_on_board = 0.0
    trapped = trapped_fraction * fuel_used
    W_oe = W_TO - fuel_used - payload - crew
    return WeightEstimate(
        W_TO=W_TO,
        fuel_used=fuel_used,
        trapped_fuel_oil=trapped,
        fuel=fuel_used + trapped,
        operating_empty_tentative=W_oe,
        empty_tentative=W_oe - trapped,
        mission_fuel_fraction=1 - fuel_used / W_TO,
    )


def roskam_empty_weight(W_TO, A, B):
    """Statistical empty weight from ``log10(W_E) = (log10(W_TO) - A) / B`` (lbf)."""
    return 10 ** ((np.log10(W_TO) - A) / B)


def size_take_off_weight(mission, payload, crew, A, B, W_guess=1e5, payload_drop_after=None, trapped_fraction=0.01):
    """Solve for the take-off weight where tentative and statistical empty weights agree."""
    from scipy.optimize import brentq

    def mismatch(W):
        est = estimate_weights(W, mission, payload, crew, payload_drop_after, trapped_fraction)
        return est.empty_tentative - roskam_empty_weight(W, A, B)

    W = brentq(mismatch, W_guess / 20, W_guess * 20)
    return estimate_weights(W, mission, payload, crew, payload_drop_after, trapped_fraction)
