"""Jackson Torok, AE 360 homework 4 (``AE360_Torok_HW4``): one day of an e = 0.1, a = 7000 km orbit
integrated with explicit Euler from periapsis.

Plots |r|, |v| and the specific energy, which drifts upward because the
Euler method does not conserve it. The MATLAB stepped 0.01 s (8.6 million
steps); the default here is 1 s (``--dt 0.01`` reproduces it, slowly).
The MATLAB also logged the energy of the state *before* each step; here
each row is consistent.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from unicodes.orbital import propagate_euler

p = argparse.ArgumentParser()
p.add_argument("--dt", type=float, default=1.0)
a_ = p.parse_args()

mu, e, a = 3.986e5, 0.1, 7000.0
rp = a * (1 - e**2) / (1 + e)
vp = np.sqrt(mu * (2 / rp - 1 / a))
t, r, v = propagate_euler([rp, 0, 0], [0, vp, 0], 86399, a_.dt, mu)
rn, vn = np.linalg.norm(r, axis=1), np.linalg.norm(v, axis=1)
E = vn**2 / 2 - mu / rn
print(f"energy: start {E[0]:.6f}, end {E[-1]:.6f} km^2/s^2 (exact {-mu / (2 * a):.6f})")
for data, label in ((rn, "Magnitude of position (km)"), (vn, "Magnitude of velocity (km/s)"),
                    (E, "Specific mechanical energy (km^2/s^2)")):
    plt.figure()
    plt.plot(t / 86400, data)
    plt.xlabel("Time (days)")
    plt.ylabel(label)
plt.show()
