"""AE 709 composites homework 4, 5a and 5b (``HW_4``, ``HW_5a``, ``HW_5b``).

AS4/3501-6 tape (psi units). HW 4 plots the transformed stiffness against
ply angle, in tension and compression moduli, with bands for +-1.5/3/4.5 deg
misalignment. HW 5a builds ABD matrices, searches the stacking permutations
of [0, +-22.5, +-45, 90] for the largest D_66 and compares E_x of three
laminates. HW 5b checks first-ply failure of [0/+-45/+-22.5/90]_s under
N_x = 500 lbf/in by maximum strain and maximum stress and finds the load
factor at first failure. (The MATLAB HW 4 computed the perturbed Q-bar at
the unperturbed angles in one loop; the perturbed angles are used.)
"""

from itertools import permutations

import matplotlib.pyplot as plt
import numpy as np

from unicodes.composites import Laminate, Ply

E_1, E_2, E_1c, G_12, nu_12 = 18.7e6, 1.9e6, 17.6e6, 0.85e6, 0.3
F = [274e3, 9.5e3, -279.6e3, -38.9e3]

# HW 4
tension, compression = Ply(E_1, E_2, G_12, nu_12), Ply(E_1c, E_2, G_12, nu_12)
angles = np.deg2rad(np.linspace(-90, 90, 181))
fig, ax = plt.subplots(figsize=(9, 5))
for (i, j), col in zip([(0, 0), (1, 1), (0, 1), (2, 2), (0, 2), (1, 2)], plt.cm.tab10.colors):
    qt = np.array([tension.Qbar(t)[i, j] for t in angles])
    qc = np.array([compression.Qbar(t)[i, j] for t in angles])
    ax.plot(np.rad2deg(angles), qt, "-", color=col, label=f"Qbar{i + 1}{j + 1} tension")
    ax.plot(np.rad2deg(angles), qc, "--", color=col)
    for pert in np.deg2rad([1.5, 3, 4.5]):
        lo = np.array([tension.Qbar(t - pert)[i, j] for t in angles])
        hi = np.array([tension.Qbar(t + pert)[i, j] for t in angles])
        ax.fill_between(np.rad2deg(angles), np.minimum(lo, hi), np.maximum(lo, hi), color=col, alpha=0.1)
ax.set(xlabel="ply angle (deg)", ylabel="psi", title="HW 4: transformed reduced stiffness")
ax.legend(fontsize=7, ncol=2)
print("HW 4: Q (tension) =\n", tension.Q)

# HW 5a
t = 0.0052
lam = Laminate([tension], None, t, np.zeros(4))
print("HW 5a [0]4: A =\n", lam.A, "\nE_x, E_y, G_xy, nu_xy =", lam.E_x, lam.E_y, lam.G_xy, lam.nu_xy)
best = max(permutations(np.deg2rad([0, 22.5, -22.5, -45, 45, 90])), key=lambda p: Laminate([tension], None, t, p).D[2, 2])
print("largest D_66 stacking:", np.round(np.rad2deg(best), 1), Laminate([tension], None, t, best).D[2, 2])
for name, ang, sym in (("[90]4", [90] * 4, 0), ("[45/-45]s", [45, -45], 2), ("[45/-45/90/90]s", [45, -45, 90, 90], 2)):
    print(f"E_x {name}: {Laminate([tension], None, t, np.deg2rad(ang), sym).E_x:.4g} psi")

# HW 5b
gamma_12 = 23000e-6  # the MATLAB wrote 23000 (micro-strain without the 1e-6)
ply = Ply(E_1, E_2, G_12, nu_12, stress_allowables=F + [G_12 * gamma_12],
          strain_allowables=[9800e-6, 4530e-6, -11300e-6, -11300e-6, gamma_12])
layup = Laminate([ply], None, t, np.deg2rad([0, 45, -45, 22.5, -22.5, 90]), 2)
loads = np.array([500.0, 0, 0, 0, 0, 0])
for criterion in ("Max_Strain", "Max_Stress"):
    failed, failures, ratios = layup.load(loads, criterion)
    r = ratios[(ratios > 0)]
    factor = r.min() if r.size else np.inf
    print(f"HW 5b {criterion}: any failure at N_x = 500: {failed.any()}; first-ply failure at N_x = {500 * factor:.0f} lbf/in")
    for f in layup.load(loads * factor * 1.0001, criterion)[1][:3]:
        print("   ", f)
plt.show()
