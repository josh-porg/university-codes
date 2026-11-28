"""AE 345 fluid mechanics: example 8.5 friction factors, homework 8 total pressures and the
closed-circuit wind-tunnel power estimate (``Wind_Tunnel_Ex8_5``, ``homework_8``, ``wind_tunnel_driver``,
``PipeSection``).
"""

import numpy as np

from unicodes.fluids import (
    PipeSection,
    darcy_pressure_drop,
    friction_factor_colebrook,
    friction_factor_haaland,
    friction_factor_laminar,
    reynolds_number,
)

# Example 8.5: 4 mm tube, 0.1 m long, air at 50 m/s
rho, mu, V, D, L = 1.23, 1.79e-5, 50.0, 4e-3, 0.1
Re = reynolds_number(rho, V, D, mu)
eod = 0.000375
print(f"Ex 8.5: Re = {Re:.0f}; laminar dp = {darcy_pressure_drop(friction_factor_laminar(Re), L, D, rho, V):.1f} Pa; "
      f"turbulent f = {friction_factor_colebrook(Re, eod):.5f} (Colebrook), {friction_factor_haaland(Re, eod):.5f} (Haaland)")
print("relative roughness  Re      Colebrook  Haaland")
for i in range(4):
    for j in range(4):
        e, R = 10.0 ** (i - 5), 10.0 ** (j + 4)
        print(f"   {e:8.0e}  {R:8.0e}  {friction_factor_colebrook(R, e):.5f}   {friction_factor_haaland(R, e):.5f}")
# (The MATLAB looped e/D = 10^(i-6) for i = 1..4 and Re = 10^(j+3); 1e3 is laminar, so Re starts at 1e4 here.)

# Homework 8: total pressure at three speeds (English units)
g, k, R, T, rho_e, p = 32.2, 1.4, 1716.0, 519.0, 2.38e-3, 2116.8
for v in (88, 330, 344.7):
    print(f"HW 8: V = {v} ft/s, M = {v / np.sqrt(k * R * T):.3f}, Bernoulli p_t = {p + 0.5 * rho_e * v**2:.1f} psf")
    # (The MATLAB added gamma*z with gamma = rho/g; z = 0 so it drops out.)

# Wind tunnel: nine sections, test section (5) at 150 ft/s, areas scaled so A_test = 6 ft^2
areas = 1.5 * np.array([22, 28, 35, 35, 4, 4, 10, 18, 22.0])
K = [0.3, 0.3, 0.3, 4.0, 0.2, 0.6, 0.3, 0.3]  # corner, corner, corner, screens, nozzle, diffuser, corner, corner
# The MATLAB indexed K from section 2 (section 1 = nan); the 8 sections between the 9 stations are used here.
V = np.empty(9)
V[4] = 150.0
for i in range(4, 8):
    V[i + 1] = PipeSection(areas[i], areas[i + 1], K[i], rho_e, V_initial=V[i]).V_final
for i in range(3, -1, -1):
    V[i] = PipeSection(areas[i], areas[i + 1], K[i], rho_e, V_final=V[i + 1]).V_initial
sections = [PipeSection(areas[i], areas[i + 1], K[i], rho_e, V_initial=V[i]) for i in range(8)]
dp = sum(s.pressure_loss for s in sections)
power = dp * V[0] * areas[0] / 0.6
print("Wind tunnel velocities (ft/s):", np.round(V, 1))
print(f"total loss {dp:.2f} psf, fan power {power / 550:.1f} hp at 60 % efficiency")
