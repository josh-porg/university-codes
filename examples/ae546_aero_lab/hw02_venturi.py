"""AE 546 Homework 2, problems 3.3-3.5: incompressible and compressible venturi.

The symbolic derivations in the MATLAB file are replaced by the closed-form
isentropic venturi result (total pressure from the two static pressures and
the area ratio).
"""

import numpy as np

from unicodes.aero import wing
from unicodes.atmosphere import isa
from unicodes.gasdynamics import isentropic

air = isa(0)

# 3.3: incompressible venturi, throat pressure
V1, area_ratio = 90.0, 1 / 0.85
dp = 0.5 * air.rho * (V1**2 - (V1 * area_ratio) ** 2)
print(f"3.3: p_throat - p_inlet = {dp:.0f} Pa (absolute {air.p + dp:.0f} Pa)")

# 3.4: mercury manometer reading 10 cm across a 12:1 venturi
g, rho_hg, dh, area_ratio = 9.81, 1.36e4, 0.1, 12.0
dp = rho_hg * g * dh
V2_incompressible = np.sqrt(2 * dp / (air.rho * (1 - (1 / area_ratio) ** 2)))
print(f"3.4 incompressible: V2 = {V2_incompressible:.1f} m/s")

gamma = 1.0005e3 / 0.718e3
p1, p2 = air.p, air.p - dp
k = (gamma + 1) / gamma
p0 = ((area_ratio**2 * p1**k - p2**k) / (area_ratio**2 * p1 ** (2 / gamma) - p2 ** (2 / gamma))) ** (
    gamma / (gamma - 1)
)
M2 = np.sqrt(2 / (gamma - 1) * ((p0 / p2) ** ((gamma - 1) / gamma) - 1))
M1 = np.sqrt(2 / (gamma - 1) * ((p0 / p1) ** ((gamma - 1) / gamma) - 1))
V2 = M2 * air.a * np.sqrt(1 / isentropic(M2, gamma).T0_T * isentropic(M1, gamma).T0_T)
print(f"3.4 compressible: p0 = {p0:.0f} Pa, M1 = {M1:.4f}, M2 = {M2:.4f}, V2 = {V2:.1f} m/s")
# The MATLAB also pushed V2 through the Karman-Tsien and Laitone formulas;
# those are pressure-coefficient corrections, so only Prandtl-Glauert is kept.
print(f"   Prandtl-Glauert correction of the incompressible V2: {wing.prandtl_glauert(V2_incompressible, M2):.1f}")

# 3.5: pitot pressure in the throat
print(f"3.5: p02 incompressible = {p2 + 0.5 * air.rho * V2_incompressible**2:.0f} Pa, isentropic = {p0:.0f} Pa")
