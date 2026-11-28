"""Jackson Torok, AE 360 homework 9 (``AE360_Torok_HW9``): Euler integration versus Kepler's solution
over one day.

From r = (1131.34, -2282.343, 6672.423) km, v = (-5.643, 4.303, 2.429) km/s
the orbit is propagated 86 399 s analytically (Kepler's equation) and with
explicit Euler steps of 0.1, 1, 10, 60 and 300 s; the table compares
positions, velocities, energy and angular momentum.

Fixes: the MATLAB's ``convert`` returned the true anomaly in degrees but
the eccentric-anomaly step used it as radians (``cos(ta)``, ``ta < pi``);
the final true anomaly added ``pi`` to a value in degrees; and the Euler
energy used ``|r|^2/2 - mu/|v|`` with r and v swapped.
"""

import numpy as np

from unicodes.orbital import propagate_euler, propagate_kepler

mu = 3.98600441e5
r0 = np.array([1131.34, -2282.343, 6672.423])
v0 = np.array([-5.64305, 4.30333, 2.4289])
T = 86399.0


def invariants(r, v):
    return np.linalg.norm(v) ** 2 / 2 - mu / np.linalg.norm(r), np.linalg.norm(np.cross(r, v))


rk, vk = propagate_kepler(r0, v0, T, mu)
Ek, Hk = invariants(rk, vk)
print(f"Kepler: r = {rk} km, v = {vk} km/s, energy {Ek:.6f}, |h| {Hk:.4f}")
for dt in (0.1, 1, 10, 60, 300):
    _, r, v = propagate_euler(r0, v0, T, dt, mu)
    E, H = invariants(r[-1], v[-1])
    print(f"Euler dt = {dt:5g} s: position error {np.linalg.norm(r[-1] - rk):10.2f} km, "
          f"velocity error {np.linalg.norm(v[-1] - vk):7.4f} km/s, energy {E:10.4f}, |h| {H:10.2f}")
