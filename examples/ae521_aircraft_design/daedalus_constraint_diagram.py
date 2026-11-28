"""Daedalus thrust-to-weight vs wing-loading constraint diagram (AE 521).

Ports ``Daedalus_preliminary_Sizing_V_0`` (AE 521 and its AE 722 copy) and
``prelim_code_to_to_run``: FAR 25 take-off and landing field lengths, the six
FAR 25 climb gradients, time to climb to 65,000 ft and cruise speed, with
Class I drag polars for each configuration. Plotted in lbf/ft^2 like the
original.
"""

import matplotlib.pyplot as plt
import numpy as np

from unicodes.aero import performance as perf
from unicodes.aero.wing import parabolic_drag_polar
from unicodes.atmosphere import RHO0, isa
from unicodes.units import FT, LBF, PSF

S_TOFL = S_FL = 8000 * FT
sigma = float(isa(2500 * FT, 57 * 5 / 9).rho) / RHO0  # 2500 ft, ISA + 57 F
W_TO = 140000 * LBF
A, e_clean, e_TO, e_L = 7.0, 0.85, 0.8, 0.75
c_f = 0.0026
W_L2W_TO = 0.4078  # from the Hermes weight sizing
WS = np.linspace(1, 500, 400) * PSF

S_wet = perf.wetted_area(W_TO, "Mil. Patrol, Bomb and Transport")
f = perf.equivalent_parasite_area(c_f, S_wet)
S = S_wet / 2
C_D0 = perf.zero_lift_drag(f, S)
print(f"S_wet = {S_wet / FT**2:.0f} ft^2, f = {f / FT**2:.1f} ft^2, C_D0 = {C_D0:.5f}")

configs = {  # name: (delta C_D0, e)
    "Clean": (0, e_clean),
    "Gear down": (0.020, e_clean),
    "Takeoff flaps": (0.01, e_TO),
    "Takeoff flaps w/ gear down": (0.03, e_TO),
    "Landing flaps": (0.055, e_L),
    "Landing flaps w/ gear down": (0.075, e_L),
}
polars = [lambda cl, d=d, e=e: parabolic_drag_polar(cl, C_D0 + d, A, e) for d, e in configs.values()]
for name, (d, e) in configs.items():
    print(f"{name:28s} C_D = {C_D0 + d:.4f} + C_L^2/{np.pi * A * e:.3f}")

fig, ax = plt.subplots()
x = WS / PSF
for cl in np.linspace(2, 1.3, 5):
    ax.plot(x, perf.t2w_takeoff_far25(S_TOFL, WS, cl, sigma), "b", lw=0.8)
for cl in np.linspace(3.8, 2.5, 5):
    ax.axvline(perf.w2s_landing_far25(S_FL, W_L2W_TO, cl, sigma) / PSF, color="g", lw=0.8)
t2w_climb = perf.t2w_climb_far25(2, WS, WS, WS, 0.975, 1.3, 2.9, polars, W_TO, 0.83 * W_TO)
for row in np.atleast_1d(t2w_climb):
    ax.axhline(row, color="r", lw=0.8)
ttc = perf.t2w_time_to_climb(3600, 65000 * FT, 70000 * FT, A, e_clean, C_D0, WS, sigma, steep=True)
ax.plot(x, ttc, "m", label="time to climb (steep)")
rho_cr = float(isa(65000 * FT).rho)
a_sl = float(isa(0).a)
ax.plot(x, perf.t2w_cruise_speed(0.5 * a_sl, A, e_clean, C_D0, WS, 1.0, rho_cr / RHO0), "k", label="cruise M 0.5 at 65 kft")
ax.set(xlabel="W/S (lbf/ft$^2$)", ylabel="T/W", xlim=(0, 500), ylim=(0, 1))
ax.legend()
plt.show()
