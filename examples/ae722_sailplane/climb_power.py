"""Climb power and high-lift sizing checks (``Climb_Power``, ``ClimbPowerCalculator``, ``HilghLift_Sizing``)."""

import numpy as np

from unicodes.aero import performance as perf
from unicodes.atmosphere import RHO0
from unicodes.units import FPM, FT, LBF, SLUG_PER_FT3

# Climb_Power (British units in the original; SI here)
W, S, rho, eta_p = 1543 * LBF, 487.3424 * FT**2, 0.002133 * SLUG_PER_FT3, 0.86
for C_D0, label in [(0.0083, "unflapped"), (0.007, "flapped")]:
    P = perf.climb_power_far(800 * FPM, eta_p, W, S, rho, C_D0, 29, 0.9)
    print(f"{label}: P = {P / 1e3:.2f} kW ({P / 2e3:.2f} kW per motor)")

# ClimbPowerCalculator: Roskam's empirical climb power
W, ws, rho = 700 * 9.81, 500.0, 1.0834
for rc_fpm in (800, 600):
    P = perf.power_for_climb(rc_fpm * FPM, eta_p, W, ws, rho / RHO0, 0.0083, 29, 0.9)
    print(f"Roskam climb power at {rc_fpm} fpm: {P / 1e3:.2f} kW")

# HilghLift_Sizing: flap C_Lmax increment needed for take-off
c_bar, taper, b, V, mu, rho = 0.85, 0.4, 20.0, 1.3 * 22.4, 1.73e-5, 1.225
c_r = 1.5 * c_bar * (taper + 1) / (taper**2 + taper + 1)
print(f"Re root {rho * V * c_r / mu:.3g}, tip {rho * V * c_r * taper / mu:.3g}")
CL_max_got = 0.95 * np.mean([1.1, 1.2]) / 1.1
dCL = 2.7 - CL_max_got
K_sweep = 1.0  # unswept
print(f"needed dC_Lmax = {dCL:.3f}, flap gives {3.4 * np.deg2rad(60) * 0.53:.3f}")
T, q = 3.13e3, 0.5 * rho * V**2
print(f"propeller slipstream dC_L = {T * np.sin(np.deg2rad(60)) / (q * c_bar * b):.3f}")
