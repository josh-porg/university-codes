"""AE 746 project 4: steady 2-D Euler flow on the project 3 grid (``Project_4_V_1``,
``configure_parameters_and_BC``, ``Initialize``, ``getQ``).

Mach-5 flow (u = 5, p = 1, rho = 1.4) over a blunt body: slip wall on the
inner boundary (unit-radius nose and afterbody), supersonic inflow through
the outer boundary, symmetry along y = 0 and supersonic outflow at the far
end. (``configure_parameters_and_BC`` listed the inflow on ``j_0`` and the
wall on ``j_n``; with the project 3 grid, where ``j = 0`` is the body,
that is the other way round.) The saved MATLAB configuration ran first
order; the second-order (minmod on conserved variables) option diverges
on this Mach-5 blunt-body case, even when restarted from the converged
first-order solution. Rusanov fluxes, SSP-RK2 with
local time steps, CFL 0.8. The grid is rebuilt as in ``project3_mesh.py``
or read from a saved ``.npz``/``.mat`` holding ``Mesh``::

    python project4_euler2d.py --order 1 --n-eta 41 --n-xi 21
    python project4_euler2d.py --mesh 2D_Mesh_641_by_321.mat --order 2 --max-iterations 60000
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from unicodes.cfd import StructuredMesh2D, euler, solve_steady
from unicodes.cfd.mesh import ruled_mesh
from unicodes.io import load_mat

p = argparse.ArgumentParser()
p.add_argument("--mesh", type=Path, default=None)
p.add_argument("--n-eta", type=int, default=41)
p.add_argument("--n-xi", type=int, default=21)
p.add_argument("--order", type=int, choices=[1, 2], default=1)
p.add_argument("--cfl", type=float, default=0.8)
p.add_argument("--global-dt", action="store_true")
p.add_argument("--max-iterations", type=int, default=10000)
a = p.parse_args()

if a.mesh is None:
    s2 = 1 / np.cos(np.pi / 6)
    s6 = (4 - s2) * np.sin(np.pi / 6)
    t = np.linspace(0, np.pi / 3, 400)
    inner = np.vstack([np.column_stack([-np.cos(t), np.sin(t)]),
                       np.linspace([-0.5, np.sqrt(3) / 2], [s6 * np.cos(np.pi / 6), s2 + s6 * np.sin(np.pi / 6)], 400)[1:]])
    t = np.linspace(0, np.pi / 6, 400)
    outer = np.vstack([np.column_stack([-2 * np.cos(t), 2 * np.sin(t)]), np.linspace([-np.sqrt(3), 1.0], [0, 4], 400)[1:]])
    nodes = ruled_mesh(inner, outer, a.n_eta, a.n_xi)
elif a.mesh.suffix == ".npz":
    nodes = np.load(a.mesh)["Mesh"]
else:
    nodes = np.asarray(load_mat(a.mesh)["Mesh"])

mesh = StructuredMesh2D(nodes)
gamma = 1.4
Q_in = euler.conserved(1.4, np.array([5.0, 0.0]), 1.0, gamma)
Q0 = np.broadcast_to(Q_in, mesh.volume.shape + (4,)).copy()
bcs = {"i_min": "symmetry", "i_max": "exit", "j_min": "wall", "j_max": ("inlet", Q_in)}
result = solve_steady(mesh, Q0, bcs, order=a.order, cfl=a.cfl, gamma=gamma, global_step=a.global_dt,
                      max_iterations=a.max_iterations)
print(f"{'converged' if result.converged else 'not converged'} after {result.iterations} iterations, "
      f"residual drop {result.residual_history[-1] / result.residual_history[0]:.2e}")

s = euler.primitives(result.Q, gamma)
X, Y = mesh.centroid[..., 0], mesh.centroid[..., 1]
fields = {"rho": s.rho, "p": s.p, "T": s.p / s.rho, "M": s.mach, "u": s.velocity[..., 0], "v": s.velocity[..., 1]}
fig, axes = plt.subplots(2, 3, figsize=(13, 8))
for ax, (name, val) in zip(axes.flat, fields.items()):
    c = ax.contourf(X, Y, val, 50, cmap="hot")
    fig.colorbar(c, ax=ax)
    ax.set(title=name, aspect="equal")
fig, ax = plt.subplots()
ax.semilogy(result.residual_history)
ax.set(title="Convergence history", xlabel="Iteration", ylabel="Density residual")
plt.show()
