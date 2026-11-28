"""AE 551 final project: longitudinal/lateral models, doublet responses and MPC pitch tracking
(``AE551_Final_Project_TF_SS_Doublet``, ``mpc1_script``, ``mpcscript_final``, ``MyPiWrap``).

The MPC Toolbox objects are replaced by :class:`unicodes.controls.LinearMPC`
with the same sample time, horizons, weights and elevator limits (+-20 deg).
The Simulink plant ``plant_C`` is not in the repository; the linear
longitudinal model stands in for it.
"""

import matplotlib.pyplot as plt
import numpy as np

from unicodes import controls as c

KT_FT = 1.68781
U1 = 112.08 * KT_FT
A_lon, B_lon = c.longitudinal_state_space(U1=1.0, theta1=np.deg2rad(6), X_u=-0.0374, X_alpha=17.4632,
                                          Z_u=-0.3404, Z_alpha=-170.4133, M_u=0, M_alpha=-5.5002, M_q=-2.1809,
                                          Z_de=-18.3007, M_de=-9.6254, Z_q=-4.0997 - 1.0)
# The project wrote the alpha row without dividing by U1 (dimensional Z terms in the A matrix); U1 = 1 and
# Z_q - 1 reproduce its matrix exactly.
A_lat, B_lat = c.lateral_state_space(U1, np.deg2rad(6), L_p=-1.0313, L_beta=-1441, L_r=0.4131, Y_p=-1.5494,
                                     Y_beta=-28.1334, Y_r=2.0192, N_p=-0.1845, N_beta=2.3097, N_r=-0.4027,
                                     L_da=2.1836, L_dr=0.2459, N_da=-0.5148, N_dr=-1.3011, Y_dr=6.02036)
for name, A in [("longitudinal", A_lon), ("lateral", A_lat)]:
    print(name, [f"{m.eigenvalue:.3f}" for m in c.modes(A)])

t = np.arange(0, 2.01, 0.01)
u = c.doublet(t, 0.5, 0.5, np.deg2rad(0.5))
fig, axes = plt.subplots(3, 4, figsize=(14, 8), sharex=True)
responses = [("elevator", c.simulate(A_lon, B_lon, u, t), ["u", "alpha", "q", "theta"]),
             ("aileron", c.simulate(A_lat, B_lat, np.column_stack([u, 0 * u]), t), ["p", "phi", "beta", "r"]),
             ("rudder", c.simulate(A_lat, B_lat, np.column_stack([0 * u, u]), t), ["p", "phi", "beta", "r"])]
for row, (name, y, labels) in zip(axes, responses):
    for ax, col, lab in zip(row, y.T, labels):
        ax.plot(t, col)
        ax.set_title(f"{lab} / {name}")

mpc = c.LinearMPC(A_lon, B_lon, C=np.array([[0, 0, 0, 1.0]]), dt=0.05, prediction_horizon=10, control_horizon=2,
                  w_y=1 * 2.7183, w_du=0.1 / 2.7183, w_u=0.0, u_min=-0.349065850398866, u_max=0.349065850398866)
X, U, Y = mpc.simulate(np.zeros(4), lambda k: np.deg2rad(5) if k > 20 else 0.0, 200)
fig, (a1, a2) = plt.subplots(2, 1, sharex=True)
k = np.arange(len(Y)) * 0.05
a1.plot(k, np.rad2deg(Y[:, 0]))
a1.set_ylabel(r"$\theta$ (deg)")
a2.step(k[1:], np.rad2deg(U[:, 0]))
a2.set_ylabel(r"$\delta_e$ (deg)")
plt.show()
