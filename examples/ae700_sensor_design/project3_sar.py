"""AE 700 project 3: airborne SAR sizing (``Project_3_algrebra``).

With the antenna length chosen so that the range and azimuth resolutions
match at the minimum look angle (``D_0 = c tau / sin(theta_min)``), solve
for the look angle at which the return from one resolution cell falls to
the minimum detectable power. ``vpasolve`` becomes ``brentq``.
"""

import numpy as np
from scipy.optimize import brentq

from unicodes import remote_sensing as rs

tau, f, P_avg, P_min = 0.01e-6, 2e9, 100.0, 1e-6
G_t = G_r = 1e3  # 30 dBi (Boeing E-7 fairing); 10**4.8 for the spacecraft option
V, H = 236.0, 1e4
lam = rs.C_LIGHT / f


def received(theta):
    D0 = rs.C_LIGHT * tau / np.sin(theta)
    R_r = rs.range_resolution(tau, theta)
    R_a = rs.azimuth_resolution(D0)
    sigma = rs.resolution_cell_rcs(R_r, R_a, theta, lam)
    rho = H / np.cos(theta)
    return rs.radar_received_power(P_avg, G_t, G_r, lam, sigma, rho, tau, rs.minimum_prf(V, D0))


theta_min = brentq(lambda t: received(t) - P_min, 0.05, 1.5)
D0 = rs.C_LIGHT * tau / np.sin(theta_min)
print(f"theta_min = {theta_min:.6f} rad ({np.degrees(theta_min):.2f} deg)")
print(f"antenna length D0 = {D0:.3f} m, PRF = {rs.minimum_prf(V, D0):.1f} Hz")
print(f"range resolution {rs.range_resolution(tau, theta_min):.3f} m, azimuth resolution {rs.azimuth_resolution(D0):.3f} m")
print(f"near-swath ground distance {H * np.tan(theta_min):.0f} m")
