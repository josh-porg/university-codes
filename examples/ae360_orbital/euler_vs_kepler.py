"""Euler integration of the two-body problem vs the analytic Kepler solution
(``Joshua_Poznanski_HW4``, ``integrate2Body``, ``driver_final``, ``ae_360_analysis``, ``Joshua_Poznanski_HW9``).

The MATLAB HW 4 used mu = 3.586e5 km^3/s^2 (Earth is 3.986e5); Earth's
value is used here. Its angular-momentum line crossed ``[x, y, z]`` with
``[z, vx, vy]``; the proper ``r x v`` is used.
"""

import numpy as np

from unicodes import orbital as o

r0 = np.array([1131.34, -2282.343, 6672.423]) * 1e3
v0 = np.array([-5.64305, 4.30333, 2.42879]) * 1e3
T = 100000.0
E0 = o.specific_energy(r0, v0)
for dt in (1, 10, 60, 300):
    t, r, v = o.propagate_euler(r0, v0, T, dt)
    dE = o.specific_energy(r[-1], v[-1]) - E0
    print(f"Euler dt = {dt:4d} s: |r| = {np.linalg.norm(r[-1]) / 1e3:10.2f} km, energy drift {dE:.4g} J/kg")
rk, vk = o.propagate_kepler(r0, v0, T)
print(f"Kepler:          |r| = {np.linalg.norm(rk) / 1e3:10.2f} km, energy drift {o.specific_energy(rk, vk) - E0:.3g} J/kg")
print("elements:", o.rv_to_coe(r0, v0))
