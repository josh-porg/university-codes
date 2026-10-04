"""Statics and dynamics group project: radii of successive roller-coaster loops with friction
(``Statics_and_dynamics_group_porject_rollercoaster``).

A 1000 kg car starts at h0 = r0 = 80 m. Friction (19.62 + 581.175 N) acts
along semicircular arcs; loop 1 is sized so the car is just weightless at
its top (v^2 = g r1), and the following hills (heights r2, r3, r4) so the car
just reaches each top. Every energy balance is linear in the unknown radius.
"""

import numpy as np

g, m, r0, fk = 9.81, 1000.0, 80.0, 19.62 + 581.175
E0 = m * g * r0
# loop 1: E0 - fk pi (r0/2 + r1) = m g r1 / 2 + m g 2 r1
r1 = (E0 - fk * np.pi * r0 / 2) / (fk * np.pi + 2.5 * m * g)
radii = [r1]
travelled = r0 / 2 + 2 * r1  # arc-length factor (divided by pi) up to the end of loop 1
for _ in range(3):
    # E0 - fk pi (travelled + r/2) = m g r
    r = (E0 - fk * np.pi * travelled) / (fk * np.pi / 2 + m * g)
    radii.append(r)
    travelled += r
for k, r in enumerate(radii, 1):
    print(f"R{k} = {r:.4f} m")
