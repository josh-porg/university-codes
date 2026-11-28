"""AE 445 Test 2, Test 3 and Quiz 4 calculations."""

import numpy as np

from unicodes.aero import wing
from unicodes.atmosphere import isa, speed_of_sound

# Test 2: pitot-static at 10 km, -43 C
T, p, pt = -43 + 274, 26436.0, 4.24e4
rho = isa(10e3).rho
V = np.sqrt(2 * (pt - p) / rho)
print(f"Test 2: V = {V:.1f} m/s, M = {V / speed_of_sound(T):.3f}")
air = isa(10e3)
q0 = 0.5 * air.rho * (0.5 * air.a) ** 2
print(f"Test 2 Q9: q = {q0:.1f} Pa, PG-corrected q = {wing.prandtl_glauert(q0, 0.5):.1f} Pa")

# Test 3 problem 4: finite-wing lift slope
print("Test 3 P4: a =", wing.finite_wing_lift_slope(5.729577951, 10, 1))

# Test 3 problem 5: model to full scale
b_m, c_m, e = 48 / 12, 6 / 12, 0.87
a_m, alpha_0 = 0.09 * 180 / np.pi, np.deg2rad(-6)
b_fs, c_fs, alpha = 49.0, 7.0, np.deg2rad(5)
a_fs = wing.rescale_lift_slope(a_m, b_m / c_m * e, b_fs / c_fs * e)
C_L = a_fs * (alpha - alpha_0)
q = 0.5 * 2.377e-3 * 150**2
print(f"Test 3 P5: a = {a_fs:.4f}/rad, C_L = {C_L:.4f}, L = {q * C_L * b_fs * c_fs:.0f} lbf")

# Test 3 problem 7: critical Mach of a 45 deg swept wing
sweep = np.deg2rad(45)
print("Test 3 P7:", wing.critical_mach_swept_2d(0.7, sweep), wing.critical_mach_swept_3d(0.7, sweep))

# Quiz 4
q0 = 0.5 * 1.225 * (0.7 * 340.3) ** 2
print("Quiz 4 P3: q =", wing.prandtl_glauert(q0, 0.7))
# S = 6000 ft^2, A = 9. (The MATLAB divided by the speed of sound instead of S.)
rho, V, W, S, a = 2.337e-3, 146.7, 600000.0, 6000.0, 1116.4
print("Quiz 4 P4: span =", np.sqrt(S * 9), " C_L =", W / (0.5 * rho * V**2 * S), " M =", V / a)
print("Quiz 4 P5: c_l =", 2 * np.pi * np.deg2rad(12 + 6))
