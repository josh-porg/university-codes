"""Schrenk spanwise lift distribution of a DC-3 at 1.1 x MGTOW (``Report_4_DC_3``)."""

import matplotlib.pyplot as plt
import numpy as np

from unicodes.aero.wing import schrenk_lift_distribution

W = 11431 * 9.81
b = 29.0
scale = (b / 2) / 7.4  # m per cm on the three-view
c_r, c_t = 2.8 * scale, 0.7 * scale
y = np.linspace(0, b / 2, 1000)
plt.plot(y, schrenk_lift_distribution(y, b, c_r, c_t, 1.1 * W))
plt.xlabel("y (m)")
plt.ylabel("Lift per span (N/m)")
plt.show()
