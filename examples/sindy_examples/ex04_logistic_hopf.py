"""SINDy example 4 (``EX04a_LogisticMap``, ``EX04b_Hopf_TVRegDiff``, ``logistic``, ``hopf``):
parameterised systems, with the bifurcation parameter as an extra state.

``--case logistic``: the noisy logistic map ``x_{k+1} = r x_k (1 - x_k)`` at
ten values of r is identified as a discrete map ``(x, r) -> (x', r)``
(fifth-order library, threshold 0.5), and the bifurcation diagrams of the
noisy data and of the identified map are drawn for r in [1, 4].

``--case hopf``: the Hopf normal form ``x' = mu x - omega y - A x (x^2 + y^2)``,
``y' = omega x + mu y - A y (x^2 + y^2)`` for eight values of mu, with
derivatives from noisy states by TV-regularised differentiation; the
identified model ``(x, y, mu)`` is integrated again for every mu.
The bifurcation diagrams iterate all r at once (the MATLAB looped over 6000
values of r one by one).
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from common import pool_data, print_model, simulate, sparse_galerkin, stls
from unicodes.numerics import tv_derivative

p = argparse.ArgumentParser()
p.add_argument("--case", choices=["logistic", "hopf"], default="logistic")
p.add_argument("--seed", type=int, default=0)
a = p.parse_args()
rng = np.random.default_rng(a.seed)


def bifurcation(step, r_values, burn=1000, keep=1000):
    """Points (x, r) after transients; each r stops once it returns within 1e-3 of its first value."""
    x = np.full(r_values.shape, 0.5)
    for _ in range(burn):
        x = step(x, r_values)
    xss = x.copy()
    active = np.ones(r_values.shape, bool)
    xs, rs = [], []
    for _ in range(keep):
        x = step(x, r_values)
        xs.append(x[active])
        rs.append(r_values[active])
        active &= ~(np.abs(x - xss) < 1e-3)
        if not active.any():
            break
    return np.concatenate(xs), np.concatenate(rs)


if a.case == "logistic":
    eps, N = 0.01, 1000
    X, DX = [], []
    for r in (2.5, 2.75, 3, 3.25, 3.5, 3.75, 3.8, 3.85, 3.9, 3.95):
        xt = np.zeros(N)
        xt[0] = 0.5
        for k in range(1, N):
            xt[k] = np.clip(r * xt[k - 1] * (1 - xt[k - 1]) + eps * rng.standard_normal(), 0, 1)
        X.append(np.column_stack([xt[:-1], np.full(N - 1, r)]))
        DX.append(np.column_stack([xt[1:], np.full(N - 1, r)]))
    X, DX = np.vstack(X), np.vstack(DX)
    Theta = pool_data(X, 5)
    Xi = stls(Theta, DX, 0.5, iterations=5)
    print_model(Xi, ["x", "r"], 5)
    print(f"relative fit residual {np.linalg.norm(Theta @ Xi - DX) / np.linalg.norm(DX):.4f}")
    a1, a2 = Xi[4, 0], Xi[7, 0]  # coefficients of x r and x^2 r
    r_values = np.arange(1, 4.0001, 0.0005)

    def noisy(x, r):
        return np.clip(r * x - r * x**2 + eps * rng.standard_normal(x.shape), 0, 1)

    fig, ax = plt.subplots(1, 2, figsize=(10, 6))
    xs, rs = bifurcation(noisy, r_values)
    ax[0].plot(xs, rs, ".", ms=1, color="k")
    ax[0].plot(X[:, 0], X[:, 1], ".", ms=3, color="r")
    xs, rs = bifurcation(lambda x, r: r * a1 * x + r * a2 * x**2, r_values)
    ax[1].plot(xs, rs, ".", ms=1, color="k")
    for axk, title in zip(ax, ("noisy data", "identified map")):
        axk.set(xlim=(0, 1), ylim=(4, 1), title=title, xlabel="x", ylabel="r")
else:
    omega, A, dt, eps = 1.0, 1.0, 0.0025, 0.005
    tspan = np.arange(1, int(round(75 / dt)) + 1) * dt

    def hopf(t, y, mu):
        r2 = y[0] ** 2 + y[1] ** 2
        return np.array([mu * y[0] - omega * y[1] - A * y[0] * r2, omega * y[0] + mu * y[1] - A * y[1] * r2])

    runs = [(mu, [2, 0], 2, 1e2) for mu in (-0.15, -0.05)]
    runs += [(mu, x0, 10, 1e1) for mu in (0.05, 0.15, 0.25, 0.35, 0.45, 0.55) for x0 in ([0.01, 0], [2, 0])]
    X, DX = [], []
    for mu, x0, alpha, ep in runs:
        xt = simulate(lambda t, y: hopf(t, y, mu), x0, tspan, rtol=1e-12, atol=1e-12)
        xt = xt + eps * rng.standard_normal(xt.shape)
        d = np.column_stack([tv_derivative(xt[:, k], 5, alpha, dx=dt, ep=ep)[:-1] for k in range(2)])
        X.append(np.column_stack([xt, np.full(len(xt), mu)])[999:-500])
        DX.append(np.column_stack([d, np.zeros(len(d))])[999:-500])
    X, DX = np.vstack(X), np.vstack(DX)
    Xi = stls(pool_data(X, 5), DX, 0.85)
    print_model(Xi, ["x", "y", "u"], 5)
    model = sparse_galerkin(Xi, 5)
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(projection="3d")
    for mu, x0, _, _ in runs:
        xi = simulate(model, [*x0, mu], tspan[::4], rtol=1e-8, atol=1e-8)
        ax.plot(xi[:, 2], xi[:, 0], xi[:, 1], lw=0.6)
    ax.set(xlabel="mu", ylabel="x", zlabel="y", title="identified Hopf normal form")
plt.show()
