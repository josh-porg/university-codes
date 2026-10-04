"""DMD book chapter 2 (``computeDMD``, ``computePOD``, ``loadDATA``, ``loadIBPM``): DMD and POD of the
cylinder wake at Re = 100.

Needs ``CYLINDER_ALL.mat`` (see :mod:`cylinder_plot`); ``--ibpm-dir``
instead assembles ``VORTALL`` from the IBPM ``ibpm00010.plt`` ...
``ibpm01500.plt`` output files as ``loadDATA`` did. DMD keeps 21 modes; POD
is computed on the data augmented with mirror images (vorticity flips sign
under reflection) to enforce the symmetric/antisymmetric mode structure::

    python ch02_fluids.py CYLINDER_ALL.mat
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from cylinder_plot import load_cylinder, plot_cylinder
from unicodes.decomposition import DMD


def load_ibpm(fname, nx=199, ny=449):
    """``(X, Y, U, V, VORT)`` of one IBPM Tecplot file (six header lines, five columns)."""
    data = np.loadtxt(fname, skiprows=6, max_rows=nx * ny)
    return [data[:, k].reshape(ny, nx) for k in range(5)]


p = argparse.ArgumentParser()
p.add_argument("data", nargs="?", help="CYLINDER_ALL.mat")
p.add_argument("--ibpm-dir", type=Path)
a = p.parse_args()
if a.ibpm_dir:
    nx, ny = 199, 449
    VORTALL = np.column_stack([load_ibpm(a.ibpm_dir / f"ibpm{10 * c:05d}.plt")[4].T.ravel(order="F") for c in range(1, 151)])
elif a.data:
    VORTALL, _, nx, ny = load_cylinder(a.data)
else:
    p.error("pass CYLINDER_ALL.mat or --ibpm-dir")

dmd = DMD(VORTALL, rank=21, truncation=None)
for i in range(9, 20, 2):
    plot_cylinder(dmd.modes[:, i].real, f"DMD mode {i + 1} (real)")
    plot_cylinder(dmd.modes[:, i].imag, f"DMD mode {i + 1} (imag)")
fig, ax = plt.subplots()
th = np.linspace(0, 2 * np.pi, 101)
ax.plot(np.cos(th), np.sin(th), "k--")
ax.scatter(dmd.eigenvalues.real, dmd.eigenvalues.imag, facecolors="none", edgecolors="k")
ax.set(xlim=(-1.1, 1.1), ylim=(-1.1, 1.1), aspect="equal", title="DMD spectrum")

# POD with mirror images
flip = np.column_stack([-np.flipud(VORTALL[:, k].reshape(ny, nx).T).T.ravel() for k in range(VORTALL.shape[1])])
Y = np.hstack([VORTALL, flip])
mean = Y.mean(axis=1)
plot_cylinder(mean, "mean wake")
PSI, S, _ = np.linalg.svd(Y - mean[:, None], full_matrices=False)
plt.figure()
plt.semilogy(S / S.sum())
plt.title("POD singular values")
for k in range(4):
    plot_cylinder(PSI[:, k] / np.abs(PSI[:, k]).max() * 5, f"POD mode {k + 1}")
plt.show()
