"""Thermal-soaring performance of a sailplane (AE 722 thermal modelling).

Drag polar ``C_D = C_D0 + C_L^2 / (pi A e)``. SI units, wing loading in N/m^2.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar

G0 = 9.80665

# Horstmann thermal model: core climb rate and its fall-off with radius
# (m/s per m) beyond 60 m.
HORSTMANN_THERMALS = {
    "A1": (1.75, 0.025),  # narrow, weak
    "A2": (3.5, 0.032),  # narrow, strong
    "B1": (1.75, 0.0045),  # wide, weak
    "B2": (3.5, 0.006),  # wide, strong
}


def horstmann_thermal(r, thermal: str = "A1"):
    """Updraft velocity (m/s) at radius ``r`` (m) for a Horstmann thermal type."""
    v_core, gradient = HORSTMANN_THERMALS[thermal]
    return v_core - gradient * (np.asarray(r, dtype=float) - 60)


def drag_coefficient(C_L, C_D0, A, e):
    return C_D0 + np.asarray(C_L) ** 2 / (np.pi * A * e)


def sink_rate(C_L, C_D0, A, e, wing_loading, rho):
    """Straight-glide sink rate (m/s, positive down) at lift coefficient ``C_L``."""
    C_L = np.asarray(C_L, dtype=float)
    V = np.sqrt(2 * wing_loading / (rho * C_L))
    return V * drag_coefficient(C_L, C_D0, A, e) / C_L


def sink_rate_in_turn(C_L, C_D0, A, e, wing_loading, rho, r, g=G0):
    """Sink rate (m/s) in a steady circling turn of radius ``r`` at ``C_L``.

    NaN where the turn is impossible at that ``C_L`` (needs more lift than level flight allows).
    """
    C_L = np.asarray(C_L, dtype=float)
    q = (2 * wing_loading / (rho * np.asarray(r, dtype=float) * g * C_L)) ** 2  # sin^2(bank)
    with np.errstate(invalid="ignore"):
        return (
            drag_coefficient(C_L, C_D0, A, e)
            * C_L ** (-1.5)
            * np.sqrt(2 * wing_loading / rho)
            * np.where(q < 1, (1 - q), np.nan) ** (-0.75)
        )


def sink_rate_at_bank(C_D, C_L, bank, wing_loading, rho):
    """Sink rate (m/s) in a turn at bank angle ``bank`` (rad) for given coefficients."""
    return C_D / (C_L**1.5 * np.cos(bank) ** 1.5) * np.sqrt(2 * wing_loading / rho)


def climb_rate(updraft, sink):
    """Net climb rate in a thermal: updraft minus the circling sink rate.

    (The MATLAB ``Vclimbspeed`` added the two, ``V_T - -V_sc``.)
    """
    return np.asarray(updraft) - np.asarray(sink)


def average_cross_country_speed(C_L, climb, C_D0, A, e, wing_loading, rho):
    """MacCready average cross-country speed (m/s).

    Glide between thermals at ``C_L``, then climb at ``climb`` (m/s) to
    regain the height lost: ``V_avg = V * climb / (climb + sink)``.
    """
    C_L = np.asarray(C_L, dtype=float)
    V = np.sqrt(2 * wing_loading / (rho * C_L))
    return V * climb / (climb + sink_rate(C_L, C_D0, A, e, wing_loading, rho))


def optimal_glide_cl(climb, C_D0, A, e, wing_loading, rho, C_L_max):
    """Inter-thermal ``C_L`` that maximises average cross-country speed.

    Returns ``(C_L_opt, V_avg_max)``.
    """
    res = minimize_scalar(
        lambda cl: -average_cross_country_speed(cl, climb, C_D0, A, e, wing_loading, rho),
        bounds=(1e-3, C_L_max),
        method="bounded",
        options={"xatol": 1e-8},
    )
    return res.x, -res.fun


def best_climb_in_thermal(thermal, C_D0, A, e, wing_loading, rho, C_L_max, r_max=300.0, g=G0):
    """Best achievable climb rate (m/s) in a Horstmann thermal, circling at ``C_L_max``.

    Returns ``(climb_rate, radius)``.
    """

    def neg_climb(r):
        s = sink_rate_in_turn(C_L_max, C_D0, A, e, wing_loading, rho, r, g)
        return np.inf if np.isnan(s) else -climb_rate(horstmann_thermal(r, thermal), s)

    r_min = 2 * wing_loading / (rho * g * C_L_max) * 1.0001  # tightest possible turn
    res = minimize_scalar(neg_climb, bounds=(r_min, r_max), method="bounded")
    return -res.fun, res.x
