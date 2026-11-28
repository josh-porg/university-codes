"""V-n (manoeuvre and gust) envelopes for FAR 23 aircraft and CS-22 gliders (AE 722).

Weights in N, areas in m^2, speeds in m/s. Each function returns the
positive and negative limit load factor at the speeds ``V`` plus the
characteristic speeds, so a V-n diagram is ``plot(V, n_pos); plot(V, n_neg)``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..atmosphere import RHO0, isa
from ..units import KMH, KT, LBF, PSF

G = 9.81


def far23_limit_load_factors(W_gross, category: str = "Normal"):
    """FAR 23.337 limit manoeuvring load factors ``(n_pos, n_neg)`` (``LoadFactor_ManuverLimit_FAR23``).

    Categories: "Normal", "Utility", "Aerobatic", "Aerobatic Glider" (JAR-22).
    (The MATLAB ``LoadFactor_ManuverLimit_FAR23`` converted lbf->N instead of
    N->lbf; both MATLAB versions gave aerobatic aircraft +4.0 and a positive
    negative limit. FAR 23.337 gives +6.0 / -3.0.)
    """
    n_pos = max(2.1 + 24000 / (W_gross / LBF + 10000), 2.5)
    if category == "Utility":
        n_pos = 4.4
    n_neg = -0.4 * n_pos
    if category == "Aerobatic":
        n_pos, n_neg = 6.0, -3.0
    if category == "Aerobatic Glider":
        n_pos, n_neg = 7.0, -5.0
    return n_pos, n_neg


@dataclass
class VnEnvelope:
    V: np.ndarray
    n_pos: np.ndarray
    n_neg: np.ndarray
    speeds: dict  # V_S1, V_A, V_C, V_D, ...


def far23_maneuver_envelope(W_gross, W_TO, S, C_L_max, V, h=0.0, category="Normal", W_stall=None):
    """FAR 23 manoeuvre envelope (``LoadFactor_Manuver_FAR23``).

    ``W_stall`` is the weight used on the stall lines (the MATLAB used the
    operating empty weight with a 1.1 factor on lift); defaults to ``W_gross``.
    Cruise speed ``V_C = k_c sqrt(W/S)`` with k_c = 33 (36 above 20 psf or aerobatic), ``V_D = 1.25 V_C``.
    """
    W_stall = W_gross if W_stall is None else W_stall
    rho = float(isa(h).rho)
    n_max, n_min = far23_limit_load_factors(W_gross, category)
    V = np.asarray(V, dtype=float)
    n_stall = 1.1 * rho * V**2 * S * C_L_max / (2 * W_stall)
    V_S1 = np.sqrt(2 * W_stall / (1.1 * rho * S * C_L_max))
    ws_psf = W_TO / S / PSF
    k_c = 33 if category in ("Normal", "Utility") and ws_psf < 20 else 36
    V_C = k_c * np.sqrt(ws_psf) * KT
    speeds = dict(V_S1=V_S1, V_A=V_S1 * np.sqrt(n_max), V_C=V_C, V_D=1.25 * V_C)
    return VnEnvelope(V, np.minimum(n_max, n_stall), np.maximum(n_min, -0.5 * n_stall), speeds)


# CS 22.337 load factors: (n1 manoeuvre +, n2 dive +, n4 manoeuvre -, n3 dive -)
CS22_LIMITS = {"Utility": (5.3, 4.0, -2.65, -1.5), "Aerobatic": (7.0, 7.0, -5.0, -5.0)}


def cs22_maneuver_envelope(W_oe, W_TO, S, C_L_max, V, h, C_D0, A, e, category="Utility", C_L_max_neg=-0.8, W_gross=None):
    """CS-22 glider manoeuvre envelope (``LoadFactor_Manuver_CS22``).

    ``V_D`` from CS 22.335(f), ``V_A = V_S1 sqrt(n1)``, ``V_NE = 0.9 V_D``,
    design gust speed ``V_B`` at least ``V_A``.
    """
    W_gross = W_TO if W_gross is None else W_gross
    rho = float(isa(h).rho)
    n1, n2, n4, n3 = CS22_LIMITS[category]
    if category == "Utility":
        V_D = 18 * (W_TO / S / 10 / C_D0) ** (1 / 3) * KMH
    else:
        V_D = (3.5 * W_gross / S / 10 + 200) * KMH
    V = np.asarray(V, dtype=float)
    pos_stall = 1.1 * rho * V**2 * S * C_L_max / (2 * W_oe)
    neg_stall = 1.1 * rho * V**2 * S * C_L_max_neg / (2 * W_oe)
    V_S1 = np.sqrt(2 * W_oe / (1.1 * rho * S * C_L_max))
    V_A = V_S1 * np.sqrt(n1)
    V_neg = np.sqrt(2 * n4 * W_oe / (1.1 * rho * S * C_L_max_neg))
    pos_high = n1 + (n2 - n1) / (V_D - V_A) * (V - V_A)
    neg_high = n4 + (n3 - n4) / (V_D - V_neg) * (V - V_neg)
    V_B = max(np.sqrt(2 * W_TO / S / (rho * np.sqrt(C_D0 * np.pi * A * e))), V_A)
    speeds = dict(V_S1=V_S1, V_A=V_A, V_B=V_B, V_D=V_D, V_T=125 * KMH, V_W=110 * KMH, V_NE=0.9 * V_D)
    n_pos = np.minimum(np.where(V <= V_A, n1, pos_high), pos_stall)
    n_neg = np.maximum(np.where(V <= V_neg, n4, neg_high), neg_stall)
    return VnEnvelope(V, n_pos, n_neg, speeds)


def gust_load_factor(V, U_gust, W, S, c_bar, C_L_alpha, rho):
    """Pratt gust load factor with the CS-22 alleviation factor (``generateGustLine``).

    ``mu_g = 2 (W/S) / (g rho c C_La)``, ``k = 0.88 mu_g / (5.3 + mu_g)``-type
    alleviation as coded in the course (``0.96 mu/(H/c) / (0.475 + mu/(H/c))``).
    """
    m = W / G
    mu = 2 * m / S / (rho * c_bar * C_L_alpha)
    H_c = 12.17 + 0.191 * mu
    k = 0.96 * (mu / H_c) / (0.475 + mu / H_c)
    return 1 + k / 2 * RHO0 * U_gust * np.asarray(V, dtype=float) * C_L_alpha / (W / S)


def cs22_gust_envelope(W, S, c_bar, C_L_alpha, V, h, V_S1, V_B, V_D):
    """CS 22.333 gust envelope: +-15 m/s at V_B, +-7.5 m/s at V_D, capped at 1.25 (V/V_S1)^2 (``LoadFactor_Gust_CS22``)."""
    rho = float(isa(h).rho)
    V = np.asarray(V, dtype=float)
    n_cap = 1.25 * (V / V_S1) ** 2

    def line(U, at=None):
        return gust_load_factor(V if at is None else at, U, W, S, c_bar, C_L_alpha, rho)

    def capped(n, cap):
        return np.where(np.abs(n) < np.abs(cap), n, cap)

    pos_C = capped(line(15.0), n_cap)
    neg_C = capped(line(-15.0), -n_cap)
    nB_pos, nB_neg, nD_pos, nD_neg = line(15.0, V_B), line(-15.0, V_B), line(7.5, V_D), line(-7.5, V_D)
    pos_high = nB_pos + (nD_pos - nB_pos) / (V_D - V_B) * (V - V_B)
    neg_high = nB_neg + (nD_neg - nB_neg) / (V_D - V_B) * (V - V_B)
    return VnEnvelope(V, np.minimum(pos_high, pos_C), np.maximum(neg_high, neg_C), dict(V_S1=V_S1, V_B=V_B, V_D=V_D))

