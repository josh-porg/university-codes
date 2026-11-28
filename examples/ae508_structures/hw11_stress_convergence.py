"""AE 508 HW 11: stress gradient extrapolation and mesh convergence for square and triangular elements.

Ports ``HW_11``, ``HW_11_gradient``, ``HW_11_tri_import`` and the plotting
helpers. The FE results are not in the repository: pass the three
``rect_*_stress*.mat`` files (each holding a ``stress`` column) or the three
NASTRAN ``tri_*_lean.txt`` element-stress listings::

    python hw11_stress_convergence.py rect rect_1_stress.mat rect_2_stress_10spaced.mat rect_3_stress_100space.mat
    python hw11_stress_convergence.py tri tri_1_lean.txt tri_2_lean.txt tri_3_lean.txt

The convergence model ``y = b (n - a) / (n - c)`` is fitted with ``curve_fit``.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from unicodes.io import load_mat


def nastran_stresses(path):
    """Element ids and max principal stresses from a NASTRAN listing (``IdAndStressFromData``)."""
    ids, stress = [], []
    for line in open(path):
        if line.startswith("0") and "SUBCASE 1 " not in line:
            values = line.split()
            ids.append(int(float(values[1])))
            stress.append(float(values[7]))
    return np.array(ids), np.array(stress)


p = argparse.ArgumentParser()
p.add_argument("kind", choices=["rect", "tri"])
p.add_argument("files", nargs=3)
args = p.parse_args()

dist = 0.5 * np.arange(6)
if args.kind == "rect":
    stresses = np.column_stack([np.ravel(load_mat(f)["stress"])[:6] for f in args.files])
    fit_rows = [slice(2, 6)] * 3
    seeds = np.array([0.5, 0.05, 0.005])
    labels = ["216", "21600", "2160000"]
else:
    (_, s1), (_, s2), (_, s3) = (nastran_stresses(f) for f in args.files)
    keep = np.array([i % 4 not in (0, 3) for i in range(s1.size)])  # MATLAB kept mod(i,4) == 2, 3 (1-based)
    s1 = s1[keep][:10]
    stresses = np.column_stack([s1[[0, 1, 2, 3, 4, 4]], s2[[0, 10, 20, 31, 40, 51]], s3[[0, 50, 99, 150, 200, 250]]])
    fit_rows = [slice(1, 5), slice(2, 6), slice(2, 6)]
    seeds = np.array([0.5, 0.05, 0.01])
    labels = ["432", "43200", "1080000"]

fits = [np.polyfit(dist[r], stresses[r, i], 1) for i, r in enumerate(fit_rows)]
fitted = np.column_stack([np.polyval(c, dist) for c in fits])
fig, ax = plt.subplots()
for i in range(3):
    line, = ax.plot(dist[1:], stresses[1:, i], "x", label=f"{labels[i]} elements")
    ax.plot(dist, fitted[:, i], color=line.get_color(), label=f"{labels[i]} elements curve fit")
ax.set(xlabel="x distance from boundary condition (in)", ylabel="Max principal stress (psi)",
       title=f"{'Square' if args.kind == 'rect' else 'Triangular'} elements stress gradient in x")
ax.legend()

nodes = (18 / seeds + 1) * (3 / seeds + 1)
dof = (nodes - 3 / seeds) * 3
edge = fitted[0]
print("Degrees of freedom:", dof, "\nExtrapolated edge stress:", edge)
try:
    (a, b, c), _ = curve_fit(lambda n, a, b, c: b * (n - a) / (n - c), dof, edge, p0=[0, edge[-1], 0], maxfev=20000)
    print(f"Convergence fit y = b (n - a)/(n - c): a = {a:.4g}, b = {b:.6g} (converged value), c = {c:.4g}")
except RuntimeError:
    print("Convergence fit did not converge")
fig, ax = plt.subplots()
ax.semilogx(dof, edge, "xr-")
ax.set(xlabel="Number of degrees of freedom", ylabel="Max principal stress (psi)", title="Convergence of stress")
plt.show()
