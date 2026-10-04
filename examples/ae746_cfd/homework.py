"""AE 746 homework 1-3: operation timing, ODE integrator and finite-difference orders, Lagrange
derivative weights and an amplification-factor check.

* HW 1 (``Homework_1``, ``HW_1_func``): time elementwise array operations.
* HW 2 (``Homework_2_numericalSolver``): y' = x y^2, y(0) = 1 to x = 1 (exact 2)
  with explicit Euler, modified Euler, SSP-RK3 and RK4 at dx = 0.1, 0.05.
  The MATLAB RK4 error used the dx = 0.1 result twice.
* HW 2 problem 3.2: forward, central and 3-point one-sided differences of
  e^x at x = 2 (the MATLAB evaluated e^(2+1) instead of e^(2+dx)).
* HW 2 problem 3.7: derivatives of the cubic Lagrange basis on four
  equally spaced points (sympy).
* ``test_solving_hw_3``: |G| <= 1 for G = (1 - i s sin b) / (1 + s^2 sin^2 b).
"""

import time

import numpy as np
import sympy as sp

from unicodes.ode import INTEGRATORS

# HW 1
n = 10**7
rng = np.random.default_rng(0)
a, b = rng.random(n), rng.random(n)
ops = {"add": lambda: a + b, "sub": lambda: a - b, "mult": lambda: a * b, "div": lambda: a / b,
       "trig": lambda: np.sin(a), "ln": lambda: np.log(a)}
print("HW 1 time per element (ns): mean / median / min over 10 runs")
for name, op in ops.items():
    times = []
    for _ in range(10):
        t0 = time.perf_counter()
        op()
        times.append(time.perf_counter() - t0)
    times = np.array(times) / n * 1e9
    print(f"  {name:5s} {times.mean():.3f} / {np.median(times):.3f} / {times.min():.3f}")

# HW 2 integrators
print("\nHW 2: y' = x y^2, error at x = 1 and observed order")
for name, method in INTEGRATORS.items():
    err = [abs(2 - method(lambda x, y: x * y**2, (0, 1), dx, 1.0)[1][-1]) for dx in (0.1, 0.05)]
    print(f"  {name:15s} errors {err[0]:.3e}, {err[1]:.3e}; order {np.log2(err[0] / err[1]):.2f}")

# HW 2 problem 3.2
print("\nHW 2 problem 3.2: d/dx e^x at x = 2")
for dx in (0.1, 0.2):
    approx = {
        "forward": (np.exp(2 + dx) - np.exp(2)) / dx,
        "central": (np.exp(2 + dx) - np.exp(2 - dx)) / (2 * dx),
        "3-point": (-3 * np.exp(2) + 4 * np.exp(2 + dx) - np.exp(2 + 2 * dx)) / (2 * dx),
    }
    print(f"  dx = {dx}:", ", ".join(f"{k} error {abs(v - np.exp(2)):.3e}" for k, v in approx.items()))

# HW 2 problem 3.7
x, h = sp.symbols("x Delta_x")
nodes = [0, h, 2 * h, 3 * h]
for i, xi in enumerate(nodes):
    others = [xj for xj in nodes if xj != xi]
    Li = sp.prod([(x - xj) / (xi - xj) for xj in others])
    print(f"L{i + 1}' at the nodes:", [sp.simplify(sp.diff(Li, x).subs(x, xj)) for xj in nodes],
          f" L{i + 1}'' at x1:", sp.simplify(sp.diff(Li, x, 2).subs(x, 0)))

# HW 3 amplification factor
s, bb = np.meshgrid(np.linspace(0, np.sqrt(3), 200), np.linspace(0, 2 * np.pi, 400))
G = (1 - 1j * s * np.sin(bb)) / (1 + s**2 * np.sin(bb) ** 2)
print(f"\nHW 3: max |G| = {np.abs(G).max():.6f} (stable if <= 1)")
