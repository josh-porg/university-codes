"""AE 573 homework 2, 3, 4, 8 and 15, final exam problem 2 and exam 2 part 2: inlet recovery,
shock trains and isentropic relations (``HW_2``, ``HW_3``, ``HW_4``, ``H_8_part_1``, ``H_8_part_2``,
``Thermo_HW_15``, ``HW_15_diffusers``, ``FinalExam_2``, ``Exam_2_part_2``).
"""

import numpy as np
from scipy.optimize import brentq

from unicodes import gasdynamics as gd

# HW 2: inlet recovery, efficiency and entropy rise
p0, M0, pt2 = 10e3, 0.85, 15.88e3
pt0 = gd.total_pressure(p0, M0)
pi_d = gd.inlet_total_pressure_recovery(pt0, pt2)
print(f"HW 2: pt0 = {pt0:.1f} Pa, pi_d = {pi_d:.4f}, eta_d = {gd.inlet_adiabatic_efficiency(M0, p0, pt2):.4f}, "
      f"Delta s / R = {gd.inlet_entropy_rise(pi_d):.4f}")


def shock_train(M, turns):
    """Oblique shocks (beta, theta in degrees) followed by a terminal normal shock."""
    recovery = 1.0
    for beta, theta in turns:
        s = gd.oblique_shock(M, np.deg2rad(beta), np.deg2rad(theta))
        M, recovery = float(s.M2), recovery * float(s.p02_p01)
    n = gd.normal_shock(M)
    return float(n.M2), recovery * float(n.p02_p01)


print("HW 3: M2 = {:.4f}, pi_d = {:.4f}".format(*shock_train(2.0, [(35, 6), (42, 8)])))
n = gd.normal_shock(3.0)
print(f"HW 4: normal shock only pi_d = {float(n.p02_p01):.4f}; two ramps + normal shock M2 = "
      "{:.4f}, pi_d = {:.4f}".format(*shock_train(3.0, [(25, 8), (28, 8)])))

# H 8 part 1: Mach number for a total/static pressure ratio of 1.276 (the MATLAB stepped M by 1e-5)
M = brentq(lambda m: float(gd.isentropic(m).p0_p) - 1.276, 0.01, 2)
print(f"H 8 part 1: M = {M:.4f}, p = {1.276 * 30e3:.0f} Pa")

# H 8 part 2: compressor shaft power. The MATLAB used Tt1 before defining it (Tt1 = Tt0 here).
gamma, e_c, pi_c, mdot, R = 1.4, 0.9, 32, 20, 287
Tt0 = 228 * float(gd.isentropic(0.8).T0_T)
P = mdot * gamma * R / (gamma - 1) * (Tt0 - Tt0 * pi_c ** ((gamma - 1) / (gamma * e_c)))
print(f"H 8 part 2: Tt0 = {Tt0:.2f} K, shaft power = {P / 1e6:.3f} MW (negative = work input)")

# Thermo HW 15: A/A* = 2
print("Thermo HW 15: M =", round(gd.mach_from_area_ratio(2, False), 5), "or", round(gd.mach_from_area_ratio(2, True), 5))

# Final exam problem 2: diffuser from M 0.8 to M 0.65
p0, T0 = 30e3, 258.15
r0, r1 = gd.isentropic(0.8), gd.isentropic(0.65)
print(f"Final exam 2: p1 = {p0 * float(r0.p0_p) / float(r1.p0_p):.1f} Pa, "
      f"T1 = {T0 * float(r0.T0_T) / float(r1.T0_T):.2f} K")

# Exam 2 part 2. M sin(beta) = 1.4 sin 40 deg = 0.90 < 1, so the 40 deg "shock" of the
# original does not exist (the relations then return a total-pressure gain).
print(f"Exam 2 part 2: normal Mach component {1.4 * np.sin(np.deg2rad(40)):.3f} (< 1: no oblique shock)")
s = gd.oblique_shock(1.4, np.deg2rad(40), np.deg2rad(6))
n = gd.normal_shock(1.1)
print(f"Exam 2 part 2: M2 behind the oblique shock = {float(s.M2):.4f}, pt ratio {float(s.p02_p01):.5f}; "
      f"normal shock at M 1.1: M2 = {float(n.M2):.4f}, pt ratio {float(n.p02_p01):.5f}")
