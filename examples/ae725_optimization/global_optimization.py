"""AE 725 global-optimisation study: CSD vs continuous simulated annealing vs a grid search.

Ports ``GlobalOptimization_CSD_Test``, ``GlobalOptimization_ContinuousSimulatedAnealing_Test``,
``ConstrainedSteepestDescentTest`` and ``Test_goldenSearch``. The HW 13.1 beam
constraints are kept but the cost ``b d`` is made multimodal with a
``1 + sin^2(2 pi b / 100)`` factor. The MATLAB gradient of that cost
repeated d f/d b in both components; the exact gradient is used here.
"""

import matplotlib.pyplot as plt
import numpy as np

from unicodes.optimize import constrained_steepest_descent, golden_section_step, inexact_step_size

# Line-search sanity checks
print("Armijo step on (x - 2)^2 from x = 1, d = 1:", 1 + inexact_step_size(lambda x: (x - 2) ** 2, 1.0, 1.0))
print("Golden-section step on (x - 3)^2 from x = 0, d = 1:", golden_section_step(lambda x: (x - 3) ** 2, 0.0, 1.0))


def fun(x):
    b, d = x
    s = np.sin(2 * np.pi * b / 100)
    f = b * (s**2 + 1) * d
    dfdb = d * (s**2 + 1) + b * d * 2 * s * np.cos(2 * np.pi * b / 100) * 2 * np.pi / 100
    return f, np.array([dfdb, b * (s**2 + 1)])


def constraint_values(b, d):
    with np.errstate(divide="ignore"):
        return np.array([(2.4e8 / (b * d**2) - 10) / 10, (4.5e5 / (2 * b * d) - 2) * 10, (d - 2 * b) / 200, -b, -d])


def constraints(x):
    b, d = x
    g = constraint_values(b, d)
    with np.errstate(divide="ignore"):
        A = np.array([[-2.4e7 / (b**2 * d**2), -4.8e7 / (b * d**3)],
                      [-2.25e6 / (b**2 * d), -2.25e6 / (b * d**2)],
                      [-1 / 100, 1 / 200], [-1, 0], [0, -1]])
    return g, A


x0 = [500, 700]
runs = {
    "CSD": constrained_steepest_descent(fun, constraints, x0, max_iterations=10000, eps1=1e-6, eps2=1e-6),
    "CSA (log cooling)": constrained_steepest_descent(fun, constraints, x0, max_iterations=1000, eps1=1e-6, eps2=1e-6,
                                                      T0=10, freeze_time=900, rng=0),
    "CSA (exponential)": constrained_steepest_descent(fun, constraints, x0, max_iterations=3000, eps1=1e-6, eps2=1e-6,
                                                      T0=100, freeze_time=2000, cooling="exponential", beta=0.99, rng=0),
}
for name, r in runs.items():
    print(f"{name}: f = {r.fun:.5g} at b = {r.x[0]:.5g}, d = {r.x[1]:.5g} ({r.iterations} iterations)")

B, D = np.meshgrid(np.linspace(1, 2000, 800), np.linspace(1, 1000, 400))
F = fun((B, D))[0]
ok = np.all(constraint_values(B, D)[:3] <= 0, axis=0)
i = np.argmin(np.where(ok, F, np.inf))
print(f"Grid ('monte carlo') minimum: f = {F.flat[i]:.5g} at b = {B.flat[i]:.5g}, d = {D.flat[i]:.5g}")

fig, ax = plt.subplots()
ax.contourf(B, D, F, 50, cmap="cool")
for k, col in zip(range(3), "rgb"):
    ax.contour(B, D, constraint_values(B, D)[k], [0], colors=col)
for (name, r), m in zip(runs.items(), "o^s"):
    h = np.array(r.x_history)
    ax.plot(h[:, 0], h[:, 1], m + "-", ms=3, lw=0.8, label=name)
ax.set(xlabel="b", ylabel="d", xlim=(0, 2000), ylim=(0, 1000), title="Iteration history")
ax.legend()
plt.show()
