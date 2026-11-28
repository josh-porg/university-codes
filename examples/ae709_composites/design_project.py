"""AE 709 design project: running loads in a composite tube under bending, axial load and torsion
(``Design_Project``, ``HIB``).

An 8 in diameter, 36 in long tube of 16 plies (0.0052 in each) carries
P_y = -8000 lbf at the tip and P_x = 65 000 lbf axially. The bending shear
flow uses the thin-walled first moment of area ``Q = (gamma - sin gamma)``
``* 2/3 (r_o^3 - r_i^3) sin(gamma)/gamma`` from the MATLAB. (``HIB``, the 3-2-1
rotation matrix, is :func:`unicodes.flight_dynamics.kinematics.dcm_inertial_to_body`.)
"""

import matplotlib.pyplot as plt
import numpy as np

P_y, P_x, d_o, L = -8000.0, 65e3, 8.0, 36.0
t = 0.0052 * 16
r_o = d_o / 2
r_i, r_m = r_o - t, r_o - t / 2
I = np.pi * r_m**3 * t
A_cs = 2 * np.pi * r_m * t


def bending_shear_flow(gamma):
    Q = (gamma - np.sin(gamma)) * 2 / 3 * (r_o**3 - r_i**3) * (np.sin(gamma) / gamma)
    return P_y * Q / I


def bending_normal_load(x, y):
    # running load = bending stress x thickness (the MATLAB returned the stress)
    return P_y * (L - x) * y / I * t


def torsional_shear(T):
    return T / (2 * np.pi * r_m**2)


gamma = np.linspace(1e-6, np.pi / 2, 100)
y = np.linspace(0, r_m, 100)
print(f"axial running load N_x = {P_x / A_cs * t:.1f} lbf/in, peak bending N_x at the root = "
      f"{bending_normal_load(0, r_m):.1f} lbf/in, peak shear flow {bending_shear_flow(np.pi / 2):.1f} lbf/in")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.plot(np.rad2deg(gamma), bending_shear_flow(gamma))
ax1.set(xlabel="gamma (deg)", ylabel="N_xy (lbf/in)", title="bending shear flow")
ax2.plot(y, bending_normal_load(0, y) + P_x / A_cs * t)
ax2.set(xlabel="y (in)", ylabel="N_x (lbf/in)", title="root normal running load")
plt.show()
