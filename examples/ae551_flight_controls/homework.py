"""AE 551 homework 3, 4 and 6: first/second-order responses, state-space modes and doublet responses.

HW 1 (log-file plotting), HW 2 (a symbolic Laplace transform) and HW 7 (an
unfinished symbolic equation) have no numerical content to port.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

from unicodes import controls as c

# HW 3: step responses
fig, ax = plt.subplots()
for num, den in [([1.2], [1, 1.2]), ([4], [1, 0.8, 4]), ([4], [1, 1.4, 4]), ([3, 9], np.polymul([1, 0.2], [4, 1])),
                 ([0.8304], [1, 0.1172, 0.8304])]:
    t, y = signal.step((num, den), T=np.linspace(0, 40, 2000))
    ax.plot(t, y, label=f"{num}/{np.round(den, 3)}")
ax.axhline(0.9, ls="--", c="m")
ax.axhline(0.95, ls="--", c="r")
ax.legend(fontsize=7)

# HW 4: longitudinal model, doublets and singlets
A = np.array([[-0.0284, 6.6202, 0, -32.1890], [-0.0002, -0.9542, 0.9942, -0.012],
              [0.0013, -7.3610, -1.3895, 0.0005], [0, 0, 1, 0]])
B = np.array([[0], [-0.0748], [-15.0071], [0]])
print("HW4 eigenvalues:", np.linalg.eigvals(A))
num, den = c.transfer_functions(A, B)
print("HW4 theta/delta_e numerator", num[3], "denominator", den)
for T in (20, 500):
    t = np.arange(0, T, 0.01)
    fig, axes = plt.subplots(4, 1, sharex=True)
    for amp in (1, 2):
        y = c.simulate(A, B, np.deg2rad(amp) * c.doublet(t, 2, 0.5), t)
        for ax, col in zip(axes, y.T):
            ax.plot(t, col)
    axes[0].set_title(f"HW4 doublet responses, {T} s")

# HW 6: Roskam dimensional derivatives -> modes
A6, B6 = c.longitudinal_state_space(U1=266.59, theta1=np.deg2rad(1.1), X_u=-0.0138, X_alpha=18.2703, Z_u=-0.1466,
                                     Z_alpha=-575.2323, M_u=0, M_alpha=-21.7904, M_q=-2.8314, Z_de=-62.6410,
                                     M_de=-23.0586, Z_q=-6.1082, Z_alpha_dot=-1.8852, M_alpha_dot=-0.7578)
for m in c.modes(A6):
    print(f"HW6 longitudinal mode {m.eigenvalue:.4f}: wn {m.natural_frequency:.4f}, zeta {m.damping_ratio:.4f}")
A_lat = np.array([[-12.9725, 0, -30.25, 2.1389], [1, 0, 0, 0], [-0.0029, 0.1463, -0.1872, -0.9917], [-0.3591, 0, 9.2719, -0.2104]])
for m in c.modes(A_lat):
    print(f"HW6 lateral mode {m.eigenvalue:.4f}: wn {m.natural_frequency:.4f}, zeta {m.damping_ratio:.4f}")
plt.show()
