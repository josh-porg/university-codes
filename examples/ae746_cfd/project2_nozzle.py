"""AE 746 project 2: quasi-1D Euler flow through a converging-diverging nozzle
(``Project_2_Main_V_4_6``, ``Project_2_MUSCL``, ``Exact_Nozzle``).

    python project2_nozzle.py --case 2 --cells 100 200 400 --order 2 --limiter vanleer
    python project2_nozzle.py --case 1 --simple-bc --local-dt

Case 1 is subsonic throughout; case 2 has a normal shock in the diverging
section. Each run is plotted against the exact solution, and the residual
histories are compared at the end (the MATLAB saved an ``Euler1DResult``
object per run; ``--out`` saves an ``.npz`` instead).
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from unicodes.cfd import euler
from unicodes.cfd.nozzle import NozzleFlow, solve_nozzle

p = argparse.ArgumentParser()
p.add_argument("--case", type=int, choices=[1, 2], default=1)
p.add_argument("--cells", type=int, nargs="+", default=[100, 200])
p.add_argument("--order", type=int, choices=[1, 2], default=2)
p.add_argument("--limiter", choices=["none", "minmod", "vanleer", "barth"], default="vanleer")
p.add_argument("--cfl", type=float, default=0.25)
p.add_argument("--local-dt", action="store_true", help="local instead of global time stepping")
p.add_argument("--simple-bc", action="store_true", help="fix boundary states to the exact solution")
p.add_argument("--out", type=Path, default=None)
a = p.parse_args()

flow = NozzleFlow(a.case)
if a.case == 2:
    print(f"Shock location x = {flow.x_shock:.5f}")
histories = {}
for n in a.cells:
    r = solve_nozzle(flow, n, order=a.order, limiter=a.limiter, cfl=a.cfl, global_step=not a.local_dt,
                     characteristic=not a.simple_bc)
    ex = flow(r.x)
    s, se = euler.primitives(r.Q), euler.primitives(ex)
    err = np.sqrt(np.mean((r.Q[:, 0] - ex[:, 0]) ** 2))
    print(f"{n} cells: converged {r.converged} in {len(r.residual_history)} iterations, RMS density error {err:.3e}")
    histories[n] = r.residual_history

    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(7, 9))
    for ax, sim, exact, label in zip(axes, (s.rho, s.velocity[:, 0], s.p, s.velocity[:, 0] / s.c),
                                     (se.rho, se.velocity[:, 0], se.p, se.velocity[:, 0] / se.c),
                                     ("rho", "u", "p", "Mach")):
        ax.plot(r.x, exact, "b-", label="exact")
        ax.plot(r.x, sim, "r.", ms=3, label="numerical")
        ax.set_ylabel(label)
    axes[0].legend()
    axes[-1].set_xlabel("x")
    fig.suptitle(f"Case {a.case}, {n} cells, order {a.order}, limiter {a.limiter}")
    if a.out:
        a.out.mkdir(parents=True, exist_ok=True)
        np.savez(a.out / f"nozzle_case{a.case}_{n}.npz", x=r.x, Q=r.Q, Q_exact=ex, residual=r.residual_history)

fig, ax = plt.subplots()
for n, h in histories.items():
    ax.semilogy(h, label=f"{n} cells")
ax.set(xlabel="Iteration", ylabel="Density residual", title="Convergence history")
ax.legend()
plt.show()
