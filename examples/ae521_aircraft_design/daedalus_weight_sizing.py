"""Daedalus weight sizing for the 400 nmi sulphuric-acid dispersal mission (``Weight_Sizing_Daedalus_V_0``).

The MATLAB read Roskam's suggested fuel fractions and empty-weight
regression from CSV files that are not in the repository; the values below
are Roskam Part I Table 2.1 (military patrol, bomb, transport) and Table 2.15
(military patrol/bomb/transport jets). Check them against your CSVs.
"""

import numpy as np

from unicodes.aero import performance as perf
from unicodes.atmosphere import isa
from unicodes.units import FT, LBF, NMI

FRACTIONS = {"start": 0.990, "taxi": 0.990, "takeoff": 0.995, "climb": 0.980, "descent": 0.990, "landing": 0.992}
A_REG, B_REG = -0.2009, 1.1037

W_dispersed = 30e3 * LBF
W_PL = W_dispersed * 1.001
I_sp = 1 / (0.34 / 3600)  # s
L_D = 16.0
V_cruise = 0.85 * float(isa(65e3 * FT).a)
V_divert = 0.85 * float(isa(20e3 * FT).a)
W_crew = perf.crew_weight(4)

W_TO = 1.6 * W_PL
for _ in range(10000):
    f = FRACTIONS
    pre = f["start"] * f["taxi"] * f["takeoff"] * f["climb"]
    cruise = perf.fuel_fraction_range(400 * NMI, V_cruise, I_sp, L_D)
    post = (f["descent"] * 1.0 * (1 + f["climb"]) / 2 * perf.fuel_fraction_range(100 * NMI, V_divert, I_sp, L_D)
            * perf.fuel_fraction_endurance(45 * 60, I_sp, L_D) * (1 + f["descent"]) / 2 * f["landing"])
    W_F = (1 - pre) * W_TO
    # cruise burn: average of starting the leg with and without the acid (the MATLAB summed these wrongly)
    W_cruise_start = pre * W_TO
    W_F += (1 - cruise) * (W_cruise_start + W_cruise_start - W_dispersed) / 2
    W_F += (1 - post) * (W_TO - W_F - W_dispersed)
    W_E_tent = W_TO - W_F - W_PL - 0.005 * W_TO - W_crew
    W_E = perf.empty_weight(W_TO, A_REG, B_REG)
    if abs(W_E_tent - W_E) / W_E < 0.005:
        break
    W_TO += (-100 if W_E_tent > W_E else 100) * LBF
print(f"W_TO = {W_TO / LBF:,.0f} lbf, W_E = {W_E_tent / LBF:,.0f} lbf, W_F = {W_F / LBF:,.0f} lbf")
print(f"mission fuel fraction pieces: pre {pre:.4f}, cruise {cruise:.4f}, post {post:.4f}")
np.set_printoptions(precision=4)
