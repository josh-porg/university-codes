"""AE 700 final exam problems 2 and 3.

Problem 2: spaceborne SAR - unresolvable region and the smallest radar
cross-section detectable at the minimum look angle. Problem 3: blackbody
spectral radiance at and around the Wien peak versus temperature. (Problem 1
is a symbolic derivation; problem 4 is a weapons-effects calculation and is
not ported.)
"""

import matplotlib.pyplot as plt
import numpy as np

from unicodes import remote_sensing as rs

# Problem 2
mu_earth, R_earth = 3.986004418e14, 6.371e6
tau, D0, f, H, L = 1e-8, 0.5, 2e9, 4e5, 4.0
P_avg, G = 100.0, 10**4.8
P_r, R_r = 1e-6, 2.0
V = np.sqrt(mu_earth / (H + R_earth))
lam = rs.C_LIGHT / f
prf = rs.minimum_prf(V, D0)
theta_min = np.arctan2(rs.C_LIGHT * tau, 2 * R_r)
rho = H * np.sin(theta_min)  # as in the MATLAB; the slant range at that look angle is H / cos(theta_min)
print(f"V = {V:.1f} m/s, synthetic aperture {H * lam / L:.1f} m, PRF = {prf:.0f} Hz")
print(f"theta_min = {theta_min:.4f} rad, unresolvable range {H * np.tan(theta_min):.0f} m")
print(f"minimum RCS {rs.rcs_for_received_power(P_r, P_avg, G, G, lam, rho, tau, prf):.4g} m^2 (rho = H sin theta), "
      f"{rs.rcs_for_received_power(P_r, P_avg, G, G, lam, H / np.cos(theta_min), tau, prf):.4g} m^2 (rho = H / cos theta)")

# Problem 3
T = np.linspace(200, 600, 100)
lam_max = rs.wien_peak_wavelength(T)
fig, ax = plt.subplots()
ax.plot(T, rs.planck_spectral_radiance(T, lam_max), "g", lw=2, label=r"$\lambda_{max}$")
ax.plot(T, rs.planck_spectral_radiance(T, 0.75 * lam_max), "b", lw=2, label=r"$0.75\lambda_{max}$")
ax.plot(T, rs.planck_spectral_radiance(T, 1.25 * lam_max), "r", lw=2, label=r"$1.25\lambda_{max}$")
ax.set(xlabel="Temperature $T$ (K)", ylabel=r"Spectral radiance $L_\lambda$ (W m$^{-2}$ sr$^{-1}$ m$^{-1}$)")
ax.legend(loc="upper left")
plt.show()
