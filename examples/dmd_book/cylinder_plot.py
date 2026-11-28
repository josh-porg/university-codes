"""Vorticity plot of the cylinder wake used by the DMD book examples (``plotCylinder``,
``plotCylinderNoSave``).

``CYLINDER_ALL.mat`` (the book's ``DATA/FLUIDS`` folder, not in this
repository) holds ``VORTALL`` (89 351 x 150 snapshots of a 199 x 449 grid),
``VORTEXTRA`` and the grid sizes ``nx``, ``ny`` (``m``, ``n``).
The book's ``CCcool`` colormap is replaced by a diverging matplotlib map.
"""

import matplotlib.pyplot as plt
import numpy as np

from unicodes.io import load_mat


def load_cylinder(path):
    d = load_mat(path)
    nx, ny = int(np.squeeze(d.get("nx", 199))), int(np.squeeze(d.get("ny", 449)))
    extra = np.asarray(d.get("VORTEXTRA", np.empty((nx * ny, 0))), float)
    return np.asarray(d["VORTALL"], float), extra.reshape(nx * ny, -1), nx, ny


def plot_cylinder(vort, title=None, ax=None):
    """``vort`` of shape (199, 449) or flat (column-major, as MATLAB reshapes)."""
    vort = np.asarray(vort)
    if vort.ndim == 1:
        vort = vort.reshape(449, 199).T
    v = np.clip(vort, -5, 5)
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 2.6))
    ax.imshow(v, cmap="RdBu_r", vmin=-5, vmax=5)
    ax.contour(v, levels=np.r_[np.arange(-5.5, -0.4, 0.5), -0.25, -0.125], colors="k", linestyles=":", linewidths=0.8)
    ax.contour(v, levels=np.r_[0.125, 0.25, np.arange(0.5, 5.6, 0.5)], colors="k", linewidths=0.8)
    th = np.linspace(0, 2 * np.pi, 100)
    ax.fill(49 + 25 * np.sin(th), 99 + 25 * np.cos(th), color=(0.3, 0.3, 0.3))
    ax.set_xticks([1, 50, 100, 150, 200, 250, 300, 350, 400, 449], ["-1", "0", "1", "2", "3", "4", "5", "6", "7", "8"])
    ax.set_yticks([1, 50, 100, 150, 199], ["2", "1", "0", "-1", "-2"])
    if title:
        ax.set_title(title)
    return ax
