"""SINDy examples 1a-1c (``EX01a_Linear2D``, ``EX01b_Cubic2D``, ``EX01c_Linear3D``): damped
oscillators identified from noisy derivatives.

1a: ``x' = A x`` with A = [[-0.1, 2], [-2, -0.1]], fifth-order library,
derivative noise 0.05, threshold 0.05; 1b: ``x' = A x^3`` (element-wise
cube); 1c: a 3-D linear system with a decaying third state, second-order
library, noise 0.01, threshold 0.085. The identified models are integrated
from the same initial condition and compared with the truth::

    python ex01_linear_cubic.py --case b
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from common import pool_data, print_model, simulate, sparse_galerkin, stls

p = argparse.ArgumentParser()
p.add_argument("--case", choices=list("abc"), default="a")
p.add_argument("--seed", type=int, default=0)
a = p.parse_args()
rng = np.random.default_rng(a.seed)

if a.case == "c":
    A = np.array([[-0.1, 2, 0], [-2, -0.1, 0], [0, 0, -0.3]])
    rhs = lambda t, x: A @ x  # noqa: E731
    polyorder, eps, lam, tspan, x0, names = 2, 0.01, 0.085, np.arange(0, 50.001, 0.01), [2, 0, 1], ["x", "y", "z"]
else:
    A = np.array([[-0.1, 2], [-2, -0.1]])
    rhs = (lambda t, x: A @ x) if a.case == "a" else (lambda t, x: A @ x**3)
    polyorder, eps, lam, tspan, x0, names = 5, 0.05, 0.05, np.arange(0, 25.001, 0.01), [2, 0], ["x", "y"]

x = simulate(rhs, x0, tspan)
dx = np.array([rhs(0, xi) for xi in x]) + eps * rng.standard_normal(x.shape)
Theta = pool_data(x, polyorder)
Xi = stls(Theta, dx, lam)
print_model(Xi, names, polyorder)

xB = simulate(sparse_galerkin(Xi, polyorder), x0, tspan)
print(f"max deviation of the identified trajectory: {np.abs(xB - x).max():.3e}")
fig = plt.figure()
if a.case == "c":
    ax = fig.add_subplot(projection="3d")
    ax.plot(*x.T, "r", lw=1.5, label="True")
    ax.plot(*xB[::5].T, "k--", lw=1.2, label="Identified")
    ax.view_init(20, 49)
else:
    ax = fig.add_subplot()
    ax.plot(x[:, 0], x[:, 1], "r", lw=1.5, label="True")
    ax.plot(xB[:, 0], xB[:, 1], "k--", lw=1.2, label="Identified")
    ax.set(xlabel="x1", ylabel="x2")
ax.legend()
plt.figure()
for k in range(x.shape[1]):
    plt.plot(tspan, x[:, k], lw=1.5, label=f"True x{k + 1}")
plt.plot(tspan[::10], xB[::10], "k--", lw=1.2)
plt.xlabel("Time")
plt.ylabel("State, x_k")
plt.legend()
plt.show()
