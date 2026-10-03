"""CS-22 manoeuvre and gust V-n diagram (``Vn_Diagram_V_0/V_1``, ``..._modified_forFAR23``)."""

import matplotlib.pyplot as plt
import numpy as np

from unicodes.aero import loads
from unicodes.aero.wing import polhamus_lift_slope

W = 700 * 9.81
A, b = 28.0, 20.0
S, c_bar = b**2 / A, b / A
C_L_alpha = polhamus_lift_slope(2 * np.pi, A, np.deg2rad(-2), 0.5 / 2.1)
C_L_max = C_L_alpha * np.deg2rad(15)
V = np.linspace(0, 110, 400)

man = loads.cs22_maneuver_envelope(W, W, S, C_L_max, V, 0.0, 0.00829, A, 0.9, "Utility")
s = man.speeds
gust = loads.cs22_gust_envelope(W, S, c_bar, C_L_alpha, V, 0.0, s["V_S1"], s["V_B"], s["V_D"])
far = loads.far23_maneuver_envelope(W, W, S, C_L_max, V, category="Utility")

fig, ax = plt.subplots()
ax.fill_between(V, man.n_neg, man.n_pos, alpha=0.1, color="b", label="CS-22 manoeuvre envelope")
ax.fill_between(V, gust.n_neg, gust.n_pos, alpha=0.1, color="g", label="CS-22 gust envelope")
ax.plot(V, far.n_pos, "k:", V, far.n_neg, "k:", label="FAR 23 utility")
for key in ("V_S1", "V_A", "V_B", "V_D", "V_NE", "V_T"):
    ax.axvline(s[key], ls="--", lw=0.8)
    ax.text(s[key], ax.get_ylim()[1] * 0.9, f" {key} = {s[key]:.3g}", rotation=90, fontsize=8)
ax.set(xlabel="V (m/s)", ylabel="n (g)")
ax.legend()
plt.show()
