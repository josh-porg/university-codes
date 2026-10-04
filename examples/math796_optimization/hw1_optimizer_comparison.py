"""Math 796 HW 1: compare optimisers on the 2-D benchmark functions (``HW1_script``, ``Homework_1``).

MATLAB solvers and their stand-ins here:

* ``fminunc`` -> BFGS
* ``patternsearch`` -> Powell (derivative-free direct search)
* ``ga`` -> :func:`unicodes.optimize.genetic_minimize`, population 20 around (20, 30)
* ``particleswarm`` -> :func:`unicodes.optimize.particle_swarm_minimize`
* ``simulannealbnd`` -> ``scipy.optimize.dual_annealing``
* ``surrogateopt`` -> ``scipy.optimize.shgo`` (no surrogate optimiser in scipy), bounds [-70, 130]
* ``ga+fminunc``, ``simulannealbnd+fminunc`` hybrids, function counts summed

All start from (20, 30). As in ``Homework_1``, each function's objective
values are normalised by the best method's value (shifted when negative) and
summarised per function category.
"""

import argparse
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import dual_annealing, minimize, shgo

from unicodes.io import to_latex_table
from unicodes.optimize import genetic_minimize, particle_swarm_minimize
from unicodes.optimize.test_functions import CATEGORIES

p = argparse.ArgumentParser()
p.add_argument("--out", type=Path, default=None, help="directory for the LaTeX summary")
args = p.parse_args()

x0 = np.array([20.0, 30.0])
lo, hi = np.array([-70.0, -70.0]), np.array([130.0, 130.0])


class Counted:
    def __init__(self, f):
        self.f, self.n = f, 0

    def __call__(self, x):
        self.n += 1
        with np.errstate(all="ignore"):
            v = float(self.f(np.asarray(x, float)))
        return v if np.isfinite(v) else 1e300


def run(name, f):
    c = Counted(f)
    rng = np.random.default_rng(0)  # MATLAB "rng default"
    if name == "fminunc":
        x = minimize(c, x0, method="BFGS").x
    elif name == "patternsearch":
        x = minimize(c, x0, method="Powell").x
    elif name == "ga":
        pop_lo, pop_hi = x0 - 10, x0 + 10
        x = genetic_minimize(c, pop_lo, pop_hi, population_size=20, generations=100, sigma=2.0, n_elites=1, rng=rng).x
    elif name == "particleswarm":
        x, _ = particle_swarm_minimize(c, lo, hi, n_particles=20, max_iterations=200, rng=rng)
    elif name == "simulannealbnd":
        x = dual_annealing(c, list(zip(lo, hi)), x0=x0, seed=0, maxiter=200).x
    elif name == "surrogateopt":
        x = shgo(c, list(zip(lo, hi)), n=64, iters=2).x
    return x, float(f(x)), c.n


base = ["fminunc", "patternsearch", "ga", "particleswarm", "simulannealbnd", "surrogateopt"]
rows = []  # (category, function, method, x, f, nfev)
for category, functions in CATEGORIES.items():
    for f in functions:
        res = {}
        for m in base:
            res[m] = run(m, f)
        for m in ("ga", "simulannealbnd"):
            c = Counted(f)
            x = minimize(c, res[m][0], method="BFGS").x
            res[m + "+fminunc"] = (x, float(f(x)), res[m][2] + c.n)
        for m, (x, fx, n) in res.items():
            rows.append((category, f.__name__, m, x, fx, n))
        print(f"{f.__name__:24s}", " ".join(f"{m}={res[m][1]:.3g}" for m in res))

# Normalise per function
norm = {}
for name in {r[1] for r in rows}:
    vals = np.array([r[4] for r in rows if r[1] == name])
    if vals.min() < 0:
        vals = vals - vals.min()
    scaled = vals / abs(vals.min()) if vals.min() != 0 else vals + 1
    for r, v in zip([r for r in rows if r[1] == name], scaled):
        norm[(r[1], r[2])] = v

methods = base + ["ga+fminunc", "simulannealbnd+fminunc"]
summary = defaultdict(dict)
for category in CATEGORIES:
    for m in methods:
        sel = [r for r in rows if r[0] == category and r[2] == m]
        summary[category][m] = (np.mean([norm[(r[1], m)] for r in sel]), np.mean([r[5] for r in sel]))

fig, axes = plt.subplots(2, 1, figsize=(10, 8))
w = 0.8 / len(methods)
for i, m in enumerate(methods):
    pos = np.arange(len(CATEGORIES)) + i * w
    axes[0].bar(pos, [summary[c][m][0] for c in CATEGORIES], w, label=m)
    axes[1].bar(pos, [summary[c][m][1] for c in CATEGORIES], w)
for ax, label in zip(axes, ["mean normalised objective", "mean function evaluations"]):
    ax.set_xticks(np.arange(len(CATEGORIES)) + 0.4, list(CATEGORIES), rotation=15)
    ax.set_yscale("log")
    ax.set_ylabel(label)
axes[0].legend(fontsize=7, ncol=2)
fig.tight_layout()

if args.out:
    args.out.mkdir(parents=True, exist_ok=True)
    to_latex_table([[c, m, *summary[c][m]] for c in CATEGORIES for m in methods],
                   ["Category", "Optimizer", "Mean normalised objective", "Mean fevals"],
                   path=args.out / "optimizer_summary.tex")
plt.show()
