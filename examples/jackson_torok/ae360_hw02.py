"""Jackson Torok, AE 360 homework 2 (``AE360_Torok_HW2``): orbital elements from a state vector and back.

Question 3 converts r = (-7953.8, -4174.5, -1009.0) km, v = (3.646,
-4.912, -4.919) km/s to classical elements; question 4 rebuilds r and v
from them (with the true anomaly rounded to a whole degree, as the
homework did) through the perifocal frame.
"""

import numpy as np

from unicodes.orbital import OrbitalElements, coe_to_rv, rv_to_coe

mu = 3.98600441e5  # km^3/s^2
r = np.array([-7953.8073, -4174.5370, -1008.9496])
v = np.array([3.6460035, -4.9118820, -4.9193608])

coe = rv_to_coe(r, v, mu)
e_vec = ((v @ v - mu / np.linalg.norm(r)) * r - (r @ v) * v) / mu
print("Question 3:")
print(f"  e vector = {e_vec}, |e| = {coe.e:.6f}")
print(f"  a = {coe.a:.4f} km, h = {np.cross(r, v)}")
print(f"  i = {np.degrees(coe.i):.4f} deg, RAAN = {np.degrees(coe.raan):.4f} deg, "
      f"argument of periapsis = {np.degrees(coe.argp):.4f} deg, true anomaly = {np.degrees(coe.nu):.4f} deg")

print("Question 4:")
nu = np.radians(round(np.degrees(coe.nu)))
R, V = coe_to_rv(OrbitalElements(coe.a, coe.e, coe.i, coe.raan, coe.argp, nu), mu)
print(f"  p = {coe.a * (1 - coe.e**2):.4f} km, true anomaly {np.degrees(nu):.0f} deg")
print(f"  R = {R} km\n  V = {V} km/s")
