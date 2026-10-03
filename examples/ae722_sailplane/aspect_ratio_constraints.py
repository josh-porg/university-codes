"""Aspect ratio vs wing loading constraints for a 20 m span (AE 722).

Sink rate in a 250 ft turn (``SinkRate_Requirment``), L/D at speed
(``Lift_to_Drag_Requirements``/``_V2``), the span limit and the combined
plot (``Combined_CrossCountry_Sizing``); also the climb-index maps of
``Climb_Coefficients``. C_D0 is scaled from the DG-1000 wetted areas.
"""

import matplotlib.pyplot as plt
import numpy as np

from unicodes.aero import performance as perf
from unicodes.aero.soaring import wing_loading_for_circling_sink
from unicodes.units import FPM, FT, KT

e, rho, b, C_L_max = 0.9, 1.0834, 20.0, 2.5
W_MGTO = 1000 * 9.81
A = np.linspace(1, 50, 100)

# C_D0 as a function of aspect ratio, scaled from the DG-1000
S_wet_fuse, S_wet_wing, S_ref, c_f = 19.4, 34.5, 17.53, 0.00245
S = b**2 / A
C_D0 = perf.zero_lift_drag(perf.equivalent_parasite_area(c_f, S_wet_fuse + S_wet_wing * S / S_ref), S) + 0.008

fig, ax = plt.subplots()
r = 250 * FT
for V_sc, colour, alpha, name in [(250 * FPM, "r", 0.3, "threshold"), (175 * FPM, "y", 0.2, "objective")]:
    ws = wing_loading_for_circling_sink(V_sc, r, C_L_max, C_D0, A, e, rho)
    ax.fill_betweenx(A, ws, 1000, color=colour, alpha=alpha, label=f"sink rate {name}")

for L_D, V, colour, alpha in [(25, 70 * KT, "y", 0.2), (30, 50 * KT, "y", 0.2), (25, 45 * KT, "r", 0.3)]:
    lo, hi = perf.w2s_for_lift_to_drag(L_D, V, C_D0, C_L_max, A, e, rho)
    lo, hi = np.nan_to_num(lo, nan=1000), np.nan_to_num(hi, nan=1000)  # unattainable: exclude the whole row
    ax.fill_betweenx(A, 0, lo, color=colour, alpha=alpha)
    ax.fill_betweenx(A, hi, 1000, color=colour, alpha=alpha, label=f"L/D {L_D} at {V / KT:.0f} kt")
    ax.axvline(C_L_max * rho * V**2 / 2, color=colour)

ax.fill_betweenx(A, 0, W_MGTO / S, color="r", alpha=0.3, label="20 m span")
ax.set(xlabel="Wing loading W/S (N/m$^2$)", ylabel="Aspect ratio A", xlim=(0, 900), ylim=(0, 50))
ax.legend(fontsize=7)

# Climb index (Roskam Part I fig. 3.23 style)
C_L_rc = perf.climb_cl_max_rate(C_D0, A, e)
index = C_L_rc**1.5 / perf.climb_cd_max_rate(C_D0)
fig, ax = plt.subplots()
sc = ax.scatter(C_L_rc, index, c=A, cmap="winter")
fig.colorbar(sc, label="A")
ax.axvline(2.5, color="k")
ax.set(xlabel=r"$C_{L_{RC,max}}$", ylabel=r"$C_L^{3/2}/C_D$")
plt.show()
