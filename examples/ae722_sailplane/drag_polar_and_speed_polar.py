"""Class I drag polars, speed polars and dive speeds for each configuration (AE 722).

From ``Class_I_Drag_Polar``, ``Speed_Polar``, ``Speed_Polars`` and
``V_sink``/``speedpolarplot``. Wetted areas use the CAD values the MATLAB
fell back on; the planform-integration and fuselage-perimeter estimates
need the ``Arychtas_*`` .mat/.csv geometry files, see
:func:`unicodes.aero.performance.wetted_area_planform` and
``wetted_area_fuselage``.
"""

import matplotlib.pyplot as plt
import numpy as np

from unicodes.aero import performance as perf
from unicodes.aero import soaring
from unicodes.aero.wing import parabolic_drag_polar
from unicodes.atmosphere import isa
from unicodes.units import FT

c_f, A, e, C_L_max = 0.00245, 29.0, 0.9, 2.5
S = 20.0**2 / A
S_wet = 15.3 + 19.8 + 6.11 + 1.38  # fuselage, wing, empennage, nacelles from CAD (m^2)
C_D0_clean = perf.zero_lift_drag(perf.equivalent_parasite_area(c_f, S_wet), S)
C_D_AB = perf.airbrake_drag_coefficient(1.5 / 7.5)
flaps, gear = 0.008, 0.02
configs = {
    "Clean": 0,
    "gear": gear,
    "Landing flaps": flaps,
    "gear & flaps L": flaps + gear,
    "airbrakes 100%": C_D_AB,
    "airbrakes 60%": 0.6 * C_D_AB,
    "airbrakes 30%": 0.3 * C_D_AB,
    "gear & 100% airbrakes": gear + C_D_AB,
    "gear, flaps & 100% airbrakes": gear + C_D_AB + flaps,
}
print(f"C_D0 clean = {C_D0_clean:.5f}")

C_L = np.linspace(0, C_L_max, 100)
rho = float(isa(3000 * FT, 10).rho)
ws = 700 * 9.81 / S
V = np.linspace(7, 70, 500)
gam = np.linspace(0, np.pi / 2, 200)
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4))
for name, dC in configs.items():
    cd0 = C_D0_clean + dC
    ax1.plot(parabolic_drag_polar(C_L, cd0, A, e), C_L, label=name)
    sink, best, minsink = soaring.speed_polar(V, cd0, A, e, ws, rho)
    ax2.plot(V, sink, label=name)
    ax2.plot(best["V"], best["sink"], "s")
    with np.errstate(invalid="ignore"):
        ax3.plot(np.rad2deg(gam), perf.dive_speed(ws, cd0, A, e, rho, gam), label=name)
ax1.set(xlabel="$C_D$", ylabel="$C_L$")
ax2.set(xlabel="V (m/s)", ylabel="Sink rate (m/s)", ylim=(-4, 0))
ax3.set(xlabel=r"Dive angle $\gamma$ (deg)", ylabel="$V_D$ (m/s)")
ax3.axvline(45, ls=":")
ax1.legend(fontsize=7)

dirtiest = C_D0_clean + max(configs.values())
_, best, _ = soaring.speed_polar(V, dirtiest, A, e, ws, rho)
if best["L_D"] < 7:
    print(f"CS-22: L/D in dirtiest configuration {best['L_D']:.3g} < 7")
print("45 deg dive speed with airbrakes:", perf.dive_speed(ws, C_D0_clean + C_D_AB, A, e, rho, np.pi / 4))
plt.show()
