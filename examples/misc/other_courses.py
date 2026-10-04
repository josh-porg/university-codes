"""Single scripts from EECS 316 and ME 712 (``EESC 316/EECS_HW_4``, ``ME 712/HW_9``).

* EECS 316 HW 4: two-node phasor circuit; the node equations
  ``v2 = (v1/z1 - i1)(z2 + z3) + v1`` and ``v1 = (v2/z4 - i2)(z2 + z3) + v2``
  are solved as a linear system for the given complex impedances and sources.
* ME 712 HW 9: Raoult's-law P-x-y diagram of an ideal binary mixture with
  saturation pressures 81.42 kPa and 39.07 kPa. (The MATLAB computed the
  pressure as ``(x_a / y_a) p_sat,a``, which equals the bubble pressure
  ``x_a p_sat,a + x_n p_sat,n``.)
"""

import matplotlib.pyplot as plt
import numpy as np

# EECS 316 HW 4
z1, z2, z3, z4 = 10, 5, 15j, -10j
i1, i2 = 1, np.exp(1j * np.deg2rad(30))
zs = z2 + z3
# v2 - (1 + zs/z1) v1 = -i1 zs ;  v1 - (1 + zs/z4) v2 = -i2 zs
M = np.array([[-(1 + zs / z1), 1], [1, -(1 + zs / z4)]])
v1, v2 = np.linalg.solve(M, [-i1 * zs, -i2 * zs])
print(f"EECS 316 HW 4: v1 = {v1:.4f} V ({abs(v1):.4f} at {np.degrees(np.angle(v1)):.2f} deg), v2 = {v2:.4f} V")

# ME 712 HW 9
p_a, p_n = 8.142e4, 3.907e4
x = np.linspace(0.001, 0.999, 1000)
p = x * p_a + (1 - x) * p_n
y = x * p_a / p
fig, ax = plt.subplots()
ax.plot(x, p / 1e3, "-b", label="bubble line (liquid)")
ax.plot(y, p / 1e3, "-r", label="dew line (vapour)")
ax.set(xlabel="mole fraction of a", ylabel="pressure (kPa)", title="ME 712 HW 9: ideal binary P-x-y")
ax.legend()
plt.show()
