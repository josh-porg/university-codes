"""Thermal modelling and average cross-country speed (AE 722 Thermal Modeling).

Horstmann thermal profiles, climb rate against turn radius for several wing
loadings (``Thermal_Modeling``), and the Quast-weighted cross-country speed
over aspect ratio and wing loading (``Themral_Modeling_Part_2``,
``untitled``, ``OptimumCL``). The MATLAB grid took ~15 min at 100x100;
``N`` is smaller here.
"""

import matplotlib.pyplot as plt
import numpy as np

from unicodes.aero import soaring

r = np.linspace(0, 200, 401)
fig, ax = plt.subplots()
for name in soaring.HORSTMANN_THERMALS:
    ax.plot(r, soaring.horstmann_thermal(r, name, core=True), label=name)
ax.set(xlabel="Thermal radius (m)", ylabel="Updraft (m/s)", ylim=(0, 6))
ax.legend()

# climb speed vs radius at C_L = 1.6, A = 20
C_D0, C_L, A, e, rho = 0.00829, 1.6, 20.0, 0.9, 1.225
radius = np.linspace(30, 175, 300)
fig, axes = plt.subplots(2, 4, figsize=(14, 6), sharex=True)
for ax, ws in zip(axes.flat, [250, 300, 350, 400, 450, 500, 550]):
    sink = soaring.sink_rate_in_turn(C_L, C_D0, A, e, ws, rho, radius)
    for name in soaring.HORSTMANN_THERMALS:
        ax.plot(radius, soaring.horstmann_thermal(radius, name) - sink, label=name)
    ax.set(title=f"W/S = {ws} N/m$^2$", ylim=(0, 4))
    ax.axvline(76.2, color="k", lw=0.5)
axes.flat[0].legend()

# Quast cross-country speed map
N = 25
e, rho, C_L_max = 0.9, 1.0834, 2.5
C_D0 = 0.0081 + 0.008  # fully flapped
wing_loadings = np.linspace(50, 900, N)
aspects = np.linspace(5, 50, N)
V = np.array([[soaring.average_speed_quast(min(1.345 * (a * e) ** 0.75 / C_D0**0.25, C_L_max), C_D0, a, e, ws, rho)
               for ws in wing_loadings] for a in aspects])
fig, ax = plt.subplots()
cs = ax.contourf(wing_loadings, aspects, V, 30, cmap="winter")
fig.colorbar(cs, label="Average cross-country speed (m/s)")
ax.set(xlabel="W/S (N/m$^2$)", ylabel="A")
plt.show()
