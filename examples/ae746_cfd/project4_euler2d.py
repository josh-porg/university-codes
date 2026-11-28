"""AE 746 project 4: steady 2-D Euler flow on the project 3 grid (``Project_4_V_1``,
``configure_parameters_and_BC``, ``Initialize``, ``getQ``).

Mach-5 flow (u = 5, p = 1, rho = 1.4) over a blunt body: slip wall on the
inner boundary (unit-radius nose and afterbody), supersonic inflow through
the outer boundary, symmetry along y = 0 and supersonic outflow at the far
end. (``configure_parameters_and_BC`` listed the inflow on ``j_0`` and the
wall on ``j_n``; with the project 3 grid, where ``j = 0`` is the body,
that is the other way round.) The saved MATLAB configuration ran first
order. Second order (``--order 2``) reconstructs the primitive variables with
minmod and falls back to first order on any face that would get a
non-positive density or pressure; limiting the conserved variables, as the
MATLAB did, diverges on this Mach-5 case. On the 41 x 21 grid the
stagnation pressure is 30.8 (first order) and 31.8 (second order) against
32.65 from the Rayleigh pitot formula. Rusanov fluxes, SSP-RK2 with
local time steps, CFL 0.8. The grid is rebuilt as in ``project3_mesh.py``
or read from a saved ``.npz``/``.mat`` holding ``Mesh``.

``--viscous`` solves the steady laminar Navier-Stokes equations instead
(no-slip wall, adiabatic unless ``--wall-temperature`` is given as a
multiple of the free-stream temperature, power-law viscosity, Pr = 0.72,
Reynolds number on the nose radius). ``--dns`` integrates the same
equations time-accurately (SSP-RK3, global time step) from an impulsive
start to ``--t-end``, with nothing modelled, and records the wall
stagnation-point pressure and temperature as the flow develops. Both
cluster the cross-wise grid points at the wall (``--stretch``) to resolve
the boundary layer; check the result by refining the grid. At Re = 1000 on
41 x 41 cells (``--stretch 2``) the steady solution converges in about
7000 iterations to a stagnation pressure of 32.9 (pitot 32.65) and an
adiabatic wall stagnation temperature of T/T_inf = 6.19 (exact 6.0); the
DNS from an impulsive start reaches 33.0 and 6.26 by t = 4 (26 000 steps,
a few minutes)::

    python project4_euler2d.py --order 1 --n-eta 41 --n-xi 21
    python project4_euler2d.py --mesh 2D_Mesh_641_by_321.mat --order 2 --max-iterations 60000
    python project4_euler2d.py --order 2 --viscous --reynolds 1000 --n-eta 41 --n-xi 41 --stretch 2
    python project4_euler2d.py --order 2 --dns --reynolds 1000 --n-eta 41 --n-xi 41 --stretch 2 --t-end 4 --save dns.npz
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from unicodes.cfd import StructuredMesh2D, Viscosity, euler, solve_steady, solve_unsteady
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
p.add_argument("--viscous", action="store_true", help="steady laminar Navier-Stokes")
p.add_argument("--dns", action="store_true", help="time-accurate Navier-Stokes (DNS)")
p.add_argument("--reynolds", type=float, default=1000.0, help="rho U R_nose / mu of the free stream")
p.add_argument("--wall-temperature", type=float, default=None, help="isothermal wall, T_w / T_inf")
p.add_argument("--stretch", type=float, default=2.0, help="wall clustering of the generated grid (NS only)")
p.add_argument("--t-end", type=float, default=4.0, help="DNS end time (nose radius / free-stream sound speed units)")
p.add_argument("--save", type=Path, default=None, help="save the DNS snapshots and monitor to this .npz")
a = p.parse_args()
viscous = a.viscous or a.dns
if a.dns and a.cfl > 0.5:
    a.cfl = 0.5  # SSP-RK3 with a global step

if a.mesh is None:
    s2 = 1 / np.cos(np.pi / 6)
    s6 = (4 - s2) * np.sin(np.pi / 6)
    t = np.linspace(0, np.pi / 3, 400)
    inner = np.vstack([np.column_stack([-np.cos(t), np.sin(t)]),
                       np.linspace([-0.5, np.sqrt(3) / 2], [s6 * np.cos(np.pi / 6), s2 + s6 * np.sin(np.pi / 6)], 400)[1:]])
    t = np.linspace(0, np.pi / 6, 400)
    outer = np.vstack([np.column_stack([-2 * np.cos(t), 2 * np.sin(t)]), np.linspace([-np.sqrt(3), 1.0], [0, 4], 400)[1:]])
    nodes = ruled_mesh(inner, outer, a.n_eta, a.n_xi, stretch=a.stretch if viscous else 0.0)
elif a.mesh.suffix == ".npz":
    nodes = np.load(a.mesh)["Mesh"]
else:
    nodes = np.asarray(load_mat(a.mesh)["Mesh"])

mesh = StructuredMesh2D(nodes)
gamma = 1.4
Q_in = euler.conserved(1.4, np.array([5.0, 0.0]), 1.0, gamma)
Q0 = np.broadcast_to(Q_in, mesh.volume.shape + (4,)).copy()
T_inf = 1.0 / 1.4  # p / (rho R) with R = 1
wall = "wall"
visc = None
if viscous:
    visc = Viscosity.from_reynolds(a.reynolds, rho=1.4, speed=5.0, length=1.0, T_ref=T_inf)
    wall = "noslip" if a.wall_temperature is None else ("noslip", a.wall_temperature * T_inf)
bcs = {"i_min": "symmetry", "i_max": "exit", "j_min": wall, "j_max": ("inlet", Q_in)}


def stagnation(Q):
    st = euler.primitives(Q[0, 0], gamma)
    return float(st.p), float(st.p / st.rho / T_inf)


if a.dns:
    out = solve_unsteady(mesh, Q0, bcs, a.t_end, order=a.order, cfl=a.cfl, gamma=gamma, viscosity=visc,
                         snapshot_every=a.t_end / 20, monitor=lambda t, Q: stagnation(Q))
    print(f"DNS to t = {out.time:.3g} in {out.steps} steps")
    Q = out.Q
    if a.save is not None:
        np.savez(a.save, nodes=nodes, Q=Q, snapshot_t=[t for t, _ in out.snapshots],
                 snapshots=np.array([q for _, q in out.snapshots]), monitor=np.array([[t, *v] for t, v in out.monitor]))
else:
    result = solve_steady(mesh, Q0, bcs, order=a.order, cfl=a.cfl, gamma=gamma, global_step=a.global_dt,
                          max_iterations=a.max_iterations, viscosity=visc)
    print(f"{'converged' if result.converged else 'not converged'} after {result.iterations} iterations, "
          f"residual drop {result.residual_history[-1] / result.residual_history[0]:.2e}")
    Q = result.Q
p_stag, T_stag = stagnation(Q)
print(f"stagnation point (first cell): p / p_inf = {p_stag:.3f} (Rayleigh pitot 32.65), T / T_inf = {T_stag:.3f}")

s = euler.primitives(Q, gamma)
X, Y = mesh.centroid[..., 0], mesh.centroid[..., 1]
fields = {"rho": s.rho, "p": s.p, "T": s.p / s.rho, "M": s.mach, "u": s.velocity[..., 0], "v": s.velocity[..., 1]}
fig, axes = plt.subplots(2, 3, figsize=(13, 8))
for ax, (name, val) in zip(axes.flat, fields.items()):
    c = ax.contourf(X, Y, val, 50, cmap="hot")
    fig.colorbar(c, ax=ax)
    ax.set(title=name, aspect="equal")
fig, ax = plt.subplots()
if a.dns:
    mon = np.array([[t, *v] for t, v in out.monitor])
    ax.plot(mon[:, 0], mon[:, 1], label="p / p_inf")
    ax.plot(mon[:, 0], mon[:, 2], label="T / T_inf")
    ax.set(title="Stagnation point history", xlabel="t")
    ax.legend()
else:
    ax.semilogy(result.residual_history)
    ax.set(title="Convergence history", xlabel="Iteration", ylabel="Density residual")
plt.show()
