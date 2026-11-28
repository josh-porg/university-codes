"""SINDy appendix A (``EXappA_Sine``): the pendulum ``x' = -sin x`` with polynomial and/or sine libraries.

Ten initial conditions in [-1.25, 1.25], dt = 0.001, t in (0, 5],
fourth-order central differences, normalised library columns. ``--case 1``
(default): fifth-order polynomials, threshold 0.1, which finds the Taylor
series of -sin x; ``--case 2``: sines only, threshold 10; ``--case 3``: both.
The identified model is integrated from the same initial conditions.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from common import central_difference, pool_data, print_model, simulate, sparse_galerkin, stls

p = argparse.ArgumentParser()
p.add_argument("--case", type=int, choices=[1, 2, 3], default=1)
a = p.parse_args()
polyorder, usesine, lam = {1: (5, False, 0.1), 2: (0, True, 10), 3: (5, True, 10)}[a.case]

dt = 0.001
tspan = np.arange(1, 5001) * dt
x0s = np.r_[np.arange(-1.25, -0.24, 0.25), np.arange(0.25, 1.26, 0.25)]


def collect(rhs):
    xs, dxs = [], []
    for x0 in x0s:
        x = simulate(rhs, [x0], tspan)[:, 0]
        xs.append(x[2:-2])
        dxs.append(central_difference(x, dt))
    return np.concatenate(xs), np.concatenate(dxs)


xall, dxall = collect(lambda t, x: -np.sin(x))
Theta = pool_data(xall[:, None], polyorder, usesine)
norms = np.linalg.norm(Theta, axis=0)
Xi = stls(Theta / norms, dxall[:, None], lam) / norms[:, None]
print_model(Xi, ["x"], polyorder, usesine)
xall2, dxall2 = collect(sparse_galerkin(Xi, polyorder, usesine))
print(f"max state error of the identified model {np.abs(xall2 - xall).max():.3e}")
fig, ax = plt.subplots(2, 1, sharex=True)
tt = dt * np.arange(1, len(xall) + 1)
ax[0].plot(tt, xall, "k", tt, xall2, "r--", lw=1.5)
ax[1].plot(tt, dxall, "k", tt, dxall2, "r--", lw=1.5)
ax[1].set_xlabel("Time")
plt.show()
