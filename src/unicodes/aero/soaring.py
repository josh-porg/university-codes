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


# Climb rate on the axis (r = 0) used for the Kubrynski parabolic core.
HORSTMANN_CORE = {"A1": 2 + 6 / 7, "A2": 5.0, "B1": 1 + 13 / 14, "B2": 3 + 6 / 7}

# Quast weather model: share of each Horstmann thermal type.
QUAST_WEIGHTS = {"A1": 0.08, "A2": 0.42, "B1": 0.08, "B2": 0.42}


def horstmann_thermal(r, thermal: str = "A1", core: bool = False):
    """Updraft velocity (m/s) at radius ``r`` (m) for a Horstmann thermal type.

    ``core=False`` is the plain linear model (``Vthermal``). ``core=True`` is
    the later ``HorstmannThermal``: inside 60 m the linear profile is averaged
    with Kubrynski's parabolic core, and the updraft is clamped at zero.
    """
    v_core, gradient = HORSTMANN_THERMALS[thermal]
    r = np.asarray(r, dtype=float)
    linear = v_core - gradient * (r - 60)
    if not core:
        return linear
    parabola = -gradient / 120 * r**2 + v_core + 30 * gradient
    return np.maximum(np.where(r > 60, linear, (linear + parabola) / 2), 0.0)


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


def wing_loading_for_circling_sink(V_sc, r, C_L, C_D0, A, e, rho, g=G0):
    """Largest wing loading (N/m^2) that circles at radius ``r`` with sink ``V_sc`` (``W2S_SinkRateCircle``).

    The MATLAB substituted into a saved symbolic solution; this solves the
    same equation numerically. Works elementwise over ``A``/``C_D0`` arrays.
    """

    def one(cd0, a):
        f = lambda ws: sink_rate_in_turn(C_L, cd0, a, e, ws, rho, r, g) - V_sc  # noqa: E731
        hi = 0.5 * rho * r * g * C_L * (1 - 1e-12)  # turn becomes impossible here
        if f(1e-6) > 0:
            return np.nan
        if f(hi * (1 - 1e-9)) < 0:
            return hi
        from scipy.optimize import brentq

        return brentq(f, 1e-6, hi * (1 - 1e-9))

    return np.vectorize(one)(C_D0, A)


def climb_in_thermal(thermal, C_L, C_D0, A, e, wing_loading, rho, r=np.arange(1, 1000, 0.1), core=True, g=G0):
    """Best climb rate and radius circling at ``C_L`` in a Horstmann thermal (grid search like the MATLAB)."""
    sink = sink_rate_in_turn(C_L, C_D0, A, e, wing_loading, rho, r, g)
    climb = horstmann_thermal(r, thermal, core) - np.where(np.isnan(sink), np.inf, sink)
    k = int(np.argmax(climb))
    return climb[k], r[k]


def interthermal_glide_cl(climb, C_D0, A, e, wing_loading, rho):
    """Optimal inter-thermal C_L from ``C_D0 - C_L^2/(pi A e) - V_c C_L^1.5 / (2 sqrt(2 W/S / rho)) = 0``.

    (``LiftCoefficient_Optimal_Interthermal_Glide``); returns 0 when there is no root.
    """
    from scipy.optimize import brentq

    f = lambda cl: C_D0 - cl**2 / (np.pi * A * e) - climb * cl**1.5 / (2 * np.sqrt(2 * wing_loading / rho))  # noqa: E731
    hi = np.sqrt(C_D0 * np.pi * A * e)
    if climb <= 0 or f(hi) > 0:
        return 0.0
    return brentq(f, 1e-9, hi)


def average_speed_horstmann(C_L_circle, C_D0, A, e, wing_loading, rho):
    """Average cross-country speed (m/s) in each Horstmann thermal (``AverageCrossCountrySpeed_Horstmnan``)."""
    out = {}
    for name in HORSTMANN_THERMALS:
        climb, _ = climb_in_thermal(name, C_L_circle, C_D0, A, e, wing_loading, rho)
        if climb <= 0:
            out[name] = 0.0
            continue
        cl = interthermal_glide_cl(climb, C_D0, A, e, wing_loading, rho)
        out[name] = float(average_cross_country_speed(cl, climb, C_D0, A, e, wing_loading, rho)) if cl > 0 else 0.0
    return out


def average_speed_quast(C_L_circle, C_D0, A, e, wing_loading, rho):
    """Quast-weather-weighted average cross-country speed (``AverageCrossCountrySpeed_Quast``)."""
    v = average_speed_horstmann(C_L_circle, C_D0, A, e, wing_loading, rho)
    return sum(QUAST_WEIGHTS[k] * v[k] for k in v)


def speed_polar(V, C_D0, A, e, wing_loading, rho):
    """Sink rate (m/s, negative down) against airspeed plus best-glide and minimum-sink points (``V_sink``).

    Returns ``(sink, best_glide, min_sink)`` where each point is a dict with
    ``V``, ``sink``, ``C_L``, ``L_D`` and flight-path angle ``gamma`` (deg).
    The MATLAB located them by finite differences; here they are analytic.
    """
    V = np.asarray(V, dtype=float)
    C_L = 2 * wing_loading / (rho * V**2)
    sink = -drag_coefficient(C_L, C_D0, A, e) / C_L * V

    def point(cl):
        v = np.sqrt(2 * wing_loading / (rho * cl))
        ld = cl / drag_coefficient(cl, C_D0, A, e)
        return dict(V=v, sink=-v / ld, C_L=cl, L_D=ld, gamma=-np.rad2deg(np.arctan(1 / ld)))

    k = np.pi * A * e
    return sink, point(np.sqrt(C_D0 * k)), point(np.sqrt(3 * C_D0 * k))


def ride_quality_index(rho, U, wing_loading, C_L_alpha, C_Y_beta, sigma_w=1.0, sigma_v=1.0):
    """Ride-discomfort index C_ride from vertical and lateral gust sensitivity (AE 722 ``Ride_quality``)."""
    a_vert = rho * U / (2 * wing_loading) * np.asarray(C_L_alpha) * sigma_w
    a_lat = abs(rho * U / (2 * wing_loading) * C_Y_beta * sigma_v)
    return np.where(a_vert > 1.6 * a_lat, 2 + 18.9 * a_vert + 12.1 * a_lat, 2 + 1.62 * a_vert + 38.9 * a_lat)
