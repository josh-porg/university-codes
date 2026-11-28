"""AE 573 exam 2 part 1 (``Exam_2_part_1``): convergent versus convergent-divergent exhaust nozzle.

Nozzle pressure ratio pt7 / p0 = 3.78, convergent-nozzle total-pressure
ratio pi_n = 0.98, Tt7 = 800 K, gamma = 1.33, cp = 1156 J/(kg K).

The MATLAB set both ``NPR = pt_7/p_0`` and ``NPR_crit = pt_7/p_0`` (so the
two could never differ), and wrote ``M_9 = (2/(gamma-1)*(pt_9/p_9)^((gamma-1)/gamma)-1)``
without the square root and with the ``- 1`` outside the ratio term.
"""

import numpy as np

from unicodes.propulsion_cycles import Gas, nozzle

gas = Gas(1.33, 1156.0)
g, R = gas.gamma, gas.R
NPR, pi_n, Tt = 3.78, 0.98, 800.0
p0 = 1.0  # everything is relative to the ambient pressure

critical = ((g + 1) / 2) ** (g / (g - 1))
print(f"critical pt9 / p9 = {critical:.4f}; critical NPR with the loss = {critical / pi_n:.4f}; NPR = {NPR}")
conv = nozzle(NPR * pi_n, Tt, p0, gas, convergent=True)
print(f"convergent: {'choked' if conv.choked else 'unchoked'}, M9 = {conv.M:.4f}, p9/p0 = {conv.p:.4f}, "
      f"T9 = {conv.T:.2f} K, V9 = {conv.V:.2f} m/s")
F_conv = conv.V + conv.pressure_thrust_per_mass(p0)
print(f"  gross thrust per unit mass flow = {F_conv:.2f} N s/kg (momentum {conv.V:.2f} + pressure "
      f"{conv.pressure_thrust_per_mass(p0):.2f})")
cd = nozzle(NPR * pi_n, Tt, p0, gas)
A_ratio = 1 / cd.M * (2 / (g + 1) * (1 + (g - 1) / 2 * cd.M**2)) ** ((g + 1) / (2 * (g - 1)))
print(f"fully expanded C-D (same loss): M9 = {cd.M:.4f}, T9 = {cd.T:.2f} K, V9 = {cd.V:.2f} m/s, A9/A8 = {A_ratio:.4f}")
print(f"  gross thrust per unit mass flow = {cd.V:.2f} N s/kg, {100 * (cd.V / F_conv - 1):.2f} % more than convergent")
print(f"speed of sound at the convergent exit = {np.sqrt(g * R * conv.T):.2f} m/s")
