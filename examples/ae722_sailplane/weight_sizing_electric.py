"""Archytas electric self-launch glider weight sizing (``Weight Sizing/Archytas_Weight_Sizing``).

Battery mass per powered segment from hard-coded Li-ion specific power
(``battPED``), iterated against a STAMPED empty-weight trend. Pass the path
of ``placemotorgliderdata.mat`` (columns: year, ?, empty mass, MTOM) to use
the trend; without it the script uses W_e/W_TO = 0.6 like the first pass.
"""

import sys

import numpy as np

from unicodes.aero import performance as perf
from unicodes.aero.sizing import stamped_trend
from unicodes.aero.wing import lift_to_drag
from unicodes.atmosphere import isa
from unicodes.units import FT

G = 9.81
eta_p, eta_elec = 0.847, 1.0
A, e, C_D0, W2S, L_D_max, V_glide = 29.0, 0.9, 0.0087, 500.0, 49.5, 29.1
M_crew, M_ballast, M_baggage, ballast_left = 245.0, 0.0, 8.0, 0.5

# Li-ion specific power (W/kg) / energy (Wh/kg) by segment duration (s), from battPED
LIION = {120: (800, 25), 15: (2200, 10), 187.5: (700, 33), 75.5: (1100, 20), 171.5: (750, 30),
         360: (500, 46), 45: (1400, 15), 562.5: (400, 54)}


def specific_power(t):
    key = min(LIION, key=lambda k: abs(k - t))
    if abs(key - t) > 1:
        raise ValueError(f"no Li-ion data for a {t:.1f} s segment (battPED only has {sorted(LIION)})")
    return LIION[key][0]


def battery_weight(W_TO, n_flights=3):
    """Mission 1: three training flights of taxi, take-off, 800 fpm climb to 2500 ft, taxi."""
    V_climb, h = 800 * FT / 60, 2500 * FT
    V = np.hypot(V_climb, 36.8)
    rho = float(isa(h / 2, 10).rho)
    W_land = W_TO - ballast_left * M_ballast * G
    segments = [  # (duration s, power W)
        (120, perf.electric_power(W_TO, 5.14, L_D_max, eta_p, eta_elec)),
        (15, perf.electric_power(W_TO, 23.15, L_D_max, eta_p, eta_elec)),
        (h / V_climb, perf.electric_power(W_TO, V, lift_to_drag(W2S, rho, V, C_D0, A, e), eta_p, eta_elec)),
        (120, perf.electric_power(W_land, 5.14, L_D_max, eta_p, eta_elec)),
    ]
    return sum(n_flights * P / specific_power(n_flights * t) for t, P in segments) * G


def battery_weight_cross_country(W_TO):
    """Mission 2: climb to 2500 ft, soar down to 1000 ft, climb to 2000 ft, then a 50 km cruise-climb."""
    V_climb = 800 * FT / 60
    V = np.hypot(V_climb, 36.8)
    W_land = W_TO - ballast_left * M_ballast * G
    gamma_glide = np.arctan(1 / L_D_max)
    gamma_climb = np.arctan(V_climb / V_glide)
    h_ccl = 50e3 / (1 / np.tan(gamma_climb) + 1 / np.tan(gamma_glide))
    t_ccl = h_ccl / np.tan(gamma_climb) / 36.8

    def climb_power(h_mid):
        return perf.electric_power(W_TO, V, lift_to_drag(W2S, float(isa(h_mid, 10).rho), V, C_D0, A, e), eta_p, eta_elec)

    segments = [
        (120, perf.electric_power(W_TO, 5.14, L_D_max, eta_p, eta_elec)),
        (15, perf.electric_power(W_TO, 23.15, L_D_max, eta_p, eta_elec)),
        (2500 * FT / V_climb, climb_power(1250 * FT)),
        (1000 * FT / V_climb, climb_power(1000 * FT)),
        (t_ccl, climb_power(1250 * FT)),
        (120, perf.electric_power(W_land, 5.14, L_D_max, eta_p, eta_elec)),
    ]
    return sum(P / specific_power(t) for t, P in segments) * G


def size(We_WTO, M_guess=831.0):
    W_payload = (M_crew + M_ballast + M_baggage) * G
    for _ in range(1000):
        W_TO = M_guess * G
        W_batt = max(battery_weight(W_TO), battery_weight_cross_country(W_TO))
        W_e_tent = W_TO - W_batt - W_payload
        W_e = We_WTO * W_TO
        if abs(W_e - W_e_tent) / W_e_tent < 0.005:
            return M_guess, W_e / G, W_batt
        M_guess += -10 if W_e < W_e_tent else 10
    raise RuntimeError("did not converge")


if __name__ == "__main__":
    ratio = 0.6
    if len(sys.argv) > 1:
        from unicodes.io import load_mat

        d = np.asarray(load_mat(sys.argv[1])["placemotorgliderdata"], dtype=float)
        ratio, _ = stamped_trend(d[:, 0], d[:, 2] / d[:, 3], 2030, n_sigma=-0.31)
    M_TO, M_e, W_batt = size(ratio)
    print(f"W_e/W_TO = {ratio:.3f}: M_TO = {M_TO:.0f} kg, M_e = {M_e:.0f} kg, battery = {W_batt / G:.1f} kg")
