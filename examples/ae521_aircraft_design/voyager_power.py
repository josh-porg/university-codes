"""Rutan Voyager power available vs required at 10,000 ft (``Report_3_voyager``)."""

import matplotlib.pyplot as plt
import numpy as np

from unicodes.aero import performance as perf
from unicodes.atmosphere import isa
from unicodes.units import FT, HP, INCH, LBF

h = 10000 * FT
rho = float(isa(h).rho)
P_fore = perf.piston_power_lapse(130 * HP, 10000 * FT, h)
P_aft = perf.piston_power_lapse(110 * HP, 10000 * FT, h)
b = 110 * FT + 8 * INCH
S = 363 * FT**2
A, e, L_D_max = b**2 / S, 0.85, 27.0
C_D0 = np.pi * A * e / (4 * L_D_max**2)
V = np.linspace(10, 120, 300)
P_req = perf.power_required(rho, V, 9695 * LBF, S, C_D0, A, e)
plt.plot(V, 0.7 * (P_fore + P_aft) * np.ones_like(V) / 1e3, label="available (eta_p = 0.7)")
plt.plot(V, P_req / 1e3, label="required at MGTOW")
plt.xlabel("V (m/s)")
plt.ylabel("Power (kW)")
plt.legend()
plt.show()
