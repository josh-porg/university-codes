"""AE 725 final project: minimum-area wing box by CSD, continuous simulated annealing, GA, PSO and SA.

Merges ``wingbox_optimization_CSD_to_submit*``, ``wingbox_optimization_continuousSA``,
``wingbox_optimization_GA``, ``ParticleSwarmOptimization(_2)``, ``Wingbox_SA_V_0``
and ``FINAL_CODE_ALL_ASSEMBLED_PART_I`` (``MasterCodeLiveScript`` is an
earlier copy of the same problem). Run e.g.::

    python wingbox_optimization.py --method csd --buckling
    python wingbox_optimization.py --method csd --stringers --epochs 20 --random-start
    python wingbox_optimization.py --method ga

The CSD iteration tables are written as LaTeX like the MATLAB
``table2latex``/``Results_Summary_Latex_Converter``.
"""

import argparse
from pathlib import Path

import numpy as np

from unicodes.io import to_latex_table
from unicodes.optimize import (
    constrained_steepest_descent,
    genetic_minimize,
    numeric_constraints,
    particle_swarm_minimize,
    simulated_annealing_minimize,
)
from unicodes.structures import WingBox

p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
p.add_argument("--method", choices=["csd", "csa", "ga", "pso", "sa"], default="csd")
p.add_argument("--buckling", action="store_true", help="add skin/web buckling constraints (bucklebox_constraints)")
p.add_argument("--stringers", action="store_true", help="make the six stringer areas design variables (stringbox)")
p.add_argument("--epochs", type=int, default=1)
p.add_argument("--random-start", action="store_true")
p.add_argument("--iterations", type=int, default=100)
p.add_argument("--T0", type=float, default=1e-3, help="continuous SA initial temperature (the MATLAB tried 1e-2 to 1e-1)")
p.add_argument("--seed", type=int, default=None)
p.add_argument("--out", type=Path, default=None, help="directory for the LaTeX iteration tables")
args = p.parse_args()

rng = np.random.default_rng(args.seed)
wb = WingBox()
n = 14 if args.stringers else 8
x_start = np.r_[1.2, 0.8, 0.8, 1.2, 0.05, 0.05, 0.05, 0.05, [0.08] * (n - 8)]
names = ["A1", "A2", "A3", "A4", "t1", "t2", "t3", "t4", "S1", "S2", "S3", "S4", "S5", "S6"][:n]


def g(x):
    return wb.constraints(x, buckling=args.buckling)


def feasible(x, tol=1e-3):
    return bool(np.all(g(x) < tol))


def penalised(x, weight=1e6):
    """Exterior penalty used by the population methods."""
    viol = g(x)
    return wb.area(x)[0] + weight * viol[viol > 0].sum()


best = None
if args.method in ("csd", "csa"):
    T0 = args.T0 if args.method == "csa" else 0.0
    freeze = args.iterations - args.iterations // 10 if args.method == "csa" else -1
    for epoch in range(args.epochs):
        x0 = rng.random(n) * 10 if args.random_start else x_start
        r = constrained_steepest_descent(wb.area, numeric_constraints(g), x0, max_iterations=args.iterations,
                                         R=1e3, T0=T0, freeze_time=freeze, rng=rng)
        print(f"Epoch {epoch + 1}: f = {r.fun:.5g}, x = {np.round(r.x, 5)}, feasible: {feasible(r.x)}")
        if feasible(r.x) and (best is None or r.fun < best.fun):
            best = r
    if best is not None and args.out:
        args.out.mkdir(parents=True, exist_ok=True)
        k = len(best.f_history)
        keep = list(range(min(5, k))) + list(range(max(5, k - 5), k))
        rows = [[i + 1, best.f_history[i], *best.x_history[i]] for i in keep]
        to_latex_table(rows, ["Iteration", "Total Area (in$^2$)", *names], path=args.out / "csd_history.tex")
        rows = [[i + 1, best.alpha_history[i], *best.d_history[i]] for i in keep if i < len(best.d_history)]
        to_latex_table(rows, ["Iteration", r"$\alpha$", *("d_" + s for s in names)], path=args.out / "csd_steps.tex")
    if best is not None:
        x_best, f_best = best.x, best.fun
elif args.method == "ga":
    r = genetic_minimize(penalised, np.zeros(n), np.r_[[2.0] * 4, [0.1] * 4, [0.2] * (n - 8)], population_size=200,
                         generations=args.iterations * 5, mutation_rate=0.05, sigma=0.01, n_elites=2, positive=True,
                         rng=rng)
    x_best, f_best = r.x, wb.area(r.x)[0]
elif args.method == "pso":
    x_best, _ = particle_swarm_minimize(penalised, np.zeros(n), np.r_[[2.0] * 4, [0.1] * 4, [0.2] * (n - 8)],
                                        n_particles=500, max_iterations=args.iterations * 3, nonnegative=True, rng=rng)
    f_best = wb.area(x_best)[0]
else:
    x_best, _ = simulated_annealing_minimize(penalised, x_start, T0=1.0, max_iterations=50_000, step=0.01,
                                             lower=1e-12, rng=rng)
    f_best = wb.area(x_best)[0]

if best is None and args.method in ("csd", "csa"):
    print("No feasible design found")
else:
    print(f"\nBest ({args.method}): area = {f_best:.5g} in^2, feasible: {feasible(x_best)}")
    for name, value in zip(names, x_best):
        print(f"  {name} = {value:.5g}")
    s = wb.stresses(x_best)
    print("  shear stresses (psi):", np.round(s.shear))
    print("  spar-cap axial stresses (psi):", np.round(s.cap_axial))
