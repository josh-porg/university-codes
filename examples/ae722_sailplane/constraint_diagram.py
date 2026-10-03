"""Archytas power-loading vs wing-loading constraint diagram (AE 722).

Combines ``Archytas_Superimposed_Design_Constraints``,
``Landing_Requirements_Reida(_edited)``, ``Takeoff_Requirements_Reworked``,
``Service_Ceiling_Requirements_Reida(_edited)`` and ``Climb_Sizing``.
ISA + 10 C throughout. (The MATLAB added 273.15 twice when building the
ISA + 10 C density; the hot-day density here comes from ``isa(h, 10)``.)
"""

import matplotlib.pyplot as plt
import numpy as np

from unicodes.aero import performance as perf
from unicodes.atmosphere import RHO0, isa
from unicodes.units import FT, FPM, LBF, PSF

WS_MAX, WP_MAX = 600.0, 0.2
ws = np.linspace(1, WS_MAX, 300)
fig, ax = plt.subplots()

# Landing: 250 m field, C_LmaxL 1.6 .. 2.4
rho_hot = float(isa(0, 10).rho)
for cl in [1.6, 1.8, 2.0, 2.2, 2.4]:
    ws_land = perf.landing_wing_loading_far23(250, cl, 1.0, rho_hot)
    ax.axvline(ws_land, label=f"Landing $C_{{L_{{max,L}}}}$ = {cl}")
    ax.axvspan(ws_land, WS_MAX, color="g", alpha=0.03)

# Service ceiling: 100 fpm at 14,000 ft, absolute ceiling 15,000 ft
L_D_max, A, e, eta_p = 46.5, 22.82, 1.0, 0.8
C_L_max = np.pi * A * e / (2 * L_D_max)
C_D0 = C_L_max**2 / (np.pi * A * e)
RC0 = perf.sea_level_rate_of_climb(100 * FPM, 14000, 15000)
sigma = float(isa(7500 * FT, 10).rho) / RHO0
wp_ceiling = perf.power_loading_climb(RC0, eta_p, 1.0, ws, sigma, C_D0, A, e)
ax.plot(ws, wp_ceiling, "r", label="Service ceiling")
ax.fill_between(ws, wp_ceiling, WP_MAX, color="r", alpha=0.1)

# Take-off: CS-22 500 m field, TOP23 = 86 lbf^2/(ft^2 hp)
sigma_to = float(isa(1000 * FT, 10).rho) / RHO0
for cl in [1.4, 1.7, 2.0]:
    wp = perf.takeoff_power_loading_far23(86, ws, cl, sigma_to)
    ax.plot(ws, wp, "--", label=f"Take-off $C_{{L_{{TO}}}}$ = {cl}")
    ax.fill_between(ws, wp, WP_MAX, color="b", alpha=0.03)

# Climb (Climb_Sizing): CS-22 minimum 90 m/min, RFP threshold 600 fpm, goal 1200 fpm
CDo = 2.44 / 480  # parasite area / wetted area estimate (ft^2/ft^2)
sigma_climb = (1.1901 + 1.1393) / 2 / RHO0
for rc, name in [(90 / 60, "CS-22 90 m/min"), (600 * FPM, "RFP 600 fpm"), (1200 * FPM, "goal 1200 fpm")]:
    ax.plot(ws, perf.power_loading_climb(rc, eta_p, 1.0, ws, sigma_climb, CDo, A, 0.9), ":", label=f"Climb {name}")

ax.axvline(2180 * LBF / (174 * FT**2), color="k", label="W/S at MTOW")
ax.plot(209, 0.1047, "p", ms=15, mfc="red", label="Design point")
ax.set(xlabel="Wing loading W/S (N/m$^2$)", ylabel="Power loading W/P (N/W)", xlim=(0, WS_MAX), ylim=(0, WP_MAX))
ax.legend(fontsize=7, ncol=3, loc="upper right")
print(f"1 psf = {PSF:.3f} Pa")
plt.show()
