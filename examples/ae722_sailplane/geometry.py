"""Wing and empennage planform layout (``wing_geometery``/``wing_geometery2``, ``Empenage_Geometery``).

Saves the corner coordinates to ``archytas_geometry.npz`` (the MATLAB saved
``Arychtas_*_Coordinates_Class_I_V_0.mat`` for ``Class_I_Drag_Polar``).
"""

import matplotlib.pyplot as plt
import numpy as np

from unicodes.aero.wing import tapered_wing, trapezoid_panel
from unicodes.units import FT, INCH

# Wing: A = 29, b = 20 m, equivalent taper 0.4, 3 m untapered inboard panel swept 5 deg (crescent layout)
A, b, taper = 29.0, 20.0, 0.4
s_i, s_o, taper_i, sweep_i = 3.0, 7.0, 1.0, np.deg2rad(5)
c_r, c_t_eq, c_bar = tapered_wing(A, b, taper)
taper_o = (taper * (s_o + s_i) - s_i * taper_i) / s_o
c_bar_i = c_bar * (s_i + s_o) / (s_i + 2 / 3 * s_o * (1 + taper_o + taper_o**2) / (1 + taper_o))
c_bar_o = 2 / 3 * c_bar_i * (1 + taper_o + taper_o**2) / (1 + taper_o)
c_t = c_bar_i * taper_o
mid = (s_i + b / 2) / 2
x_w = np.array([0, s_i, mid, b / 2, b / 2, mid, s_i, 0])
y_w = np.array([0, mid * np.sin(sweep_i), (mid * np.sin(sweep_i) + c_t) / 2 + c_t, 0, -c_t,
                mid * np.sin(sweep_i) - c_bar_o + c_t, mid * np.sin(sweep_i) - c_bar_i, -c_bar_i])
print(f"equivalent wing: c_r = {c_r:.3f} m, c_bar = {c_bar:.3f} m, outboard taper = {taper_o:.3f}")

b_h, x_h, y_h, mgc_h = trapezoid_panel(14.3 * FT**2, 16 / 28, 20 * INCH, np.deg2rad(15.75))
b_v, x_v, y_v, mgc_v = trapezoid_panel(14 * FT**2, 20 / 30, 30 * INCH, np.deg2rad(35), full_span=False)
print(f"horizontal tail span {b_h:.3f} m, MGC {mgc_h[0]:.3f} m; vertical tail height {b_v:.3f} m, MGC {mgc_v[0]:.3f} m")
np.savez("archytas_geometry.npz", wing=np.vstack([x_w, y_w]), htail=np.vstack([x_h, y_h]), htail_mgc=mgc_h,
         vtail=np.vstack([x_v, y_v]), vtail_mgc=mgc_v)

fig, ax = plt.subplots()
ax.plot(np.r_[x_w, -x_w[::-1]], np.r_[y_w, y_w[::-1]], "-x", label="wing")
ax.plot(np.r_[x_h, -x_h[::-1]], np.r_[y_h, y_h[::-1]] - 5, "-x", label="horizontal tail (offset)")
ax.set_aspect("equal")
ax.legend()
plt.show()
