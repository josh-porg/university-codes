"""AE 746 project 3: algebraic grid between two arc-plus-line boundaries (``project_3_main_v_1``,
``arcLenParamterizationPlaybox``, ``thirdArc``).

The inner boundary is a 60 deg arc of radius 1 followed by a straight line,
the outer a 30 deg arc of radius 2 followed by a line; both are distributed
by arc length and joined by straight grid lines. The MATLAB built the
arc-length parameterisation by hand (and used sqrt(3/2) for sqrt(3)/2 in the
inner length); the polylines are resampled numerically here. Saves
``2D_Mesh_<n_eta>_by_<n_xi>.npz`` with the node array.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from unicodes.cfd import StructuredMesh2D
from unicodes.cfd.mesh import ruled_mesh

p = argparse.ArgumentParser()
p.add_argument("--n-eta", type=int, default=41)
p.add_argument("--n-xi", type=int, default=21)
p.add_argument("--save", action="store_true")
a = p.parse_args()

s2 = 1 / np.cos(np.pi / 6)
s6 = (4 - s2) * np.sin(np.pi / 6)
inner_end = np.array([s6 * np.cos(np.pi / 6), s2 + s6 * np.sin(np.pi / 6)])
outer_end = np.array([0.0, 4.0])

t = np.linspace(0, np.pi / 3, 400)  # from (-1, 0) up to 60 deg
inner = np.vstack([np.column_stack([-np.cos(t), np.sin(t)]), np.linspace([-0.5, np.sqrt(3) / 2], inner_end, 400)[1:]])
t = np.linspace(0, np.pi / 6, 400)  # from (-2, 0) up to 30 deg
outer = np.vstack([np.column_stack([-2 * np.cos(t), 2 * np.sin(t)]), np.linspace([-np.sqrt(3), 1.0], outer_end, 400)[1:]])

nodes = ruled_mesh(inner, outer, a.n_eta, a.n_xi)
mesh = StructuredMesh2D(nodes)
print(f"{a.n_eta} x {a.n_xi} nodes, cell areas {mesh.volume.min():.4g} to {mesh.volume.max():.4g}")

fig, ax = plt.subplots()
ax.plot(nodes[..., 0], nodes[..., 1], "b-", lw=0.6)
ax.plot(nodes[..., 0].T, nodes[..., 1].T, "b-", lw=0.6)
ax.plot(inner[:, 0], inner[:, 1], "k", outer[:, 0], outer[:, 1], "k")
ax.set_aspect("equal")
ax.set_title("Project 3 mesh")
if a.save:
    np.savez(f"2D_Mesh_{a.n_eta}_by_{a.n_xi}.npz", Mesh=nodes)
plt.show()
