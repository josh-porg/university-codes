"""SINDy example 2 (``EX02_Lorenz``, ``EX02_LorenzTVDiff``, ``lorenz``, ``TVRegDiff``): the Lorenz system.

Default: exact derivatives with noise of strength 1 added, fifth-order
library, threshold 0.025, data for t in (0, 100] at dt = 0.001.
``--tvdiff``: noise of 0.01 on the *states*, derivatives by total-variation
regularised differentiation (:func:`unicodes.numerics.tv_derivative`), the
states re-integrated from those derivatives and the ends trimmed, t in (0, 50].
The identified model is integrated to t = 20 and t = 250 next to the truth
(the attractor matches though the chaotic trajectories separate).
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from common import color_line3, pool_data, print_model, simulate, sparse_galerkin, stls
from unicodes.decomposition import lorenz
from unicodes.numerics import tv_derivative

p = argparse.ArgumentParser()
p.add_argument("--tvdiff", action="store_true")
p.add_argument("--seed", type=int, default=0)
a = p.parse_args()
rng = np.random.default_rng(a.seed)

polyorder, x0, dt = 5, np.array([-8.0, 8, 27]), 0.001
T_end = 50 if a.tvdiff else 100
tspan = np.arange(1, int(round(T_end / dt)) + 1) * dt
x = simulate(lorenz, x0, tspan, rtol=1e-12, atol=1e-12)
dxclean = np.array([lorenz(0, xi) for xi in x])
if a.tvdiff:
    xn = x + 0.01 * rng.standard_normal(x.shape)
    dxt = np.column_stack([tv_derivative(xn[:, k], 10, 2e-5, dx=dt, ep=1e12)[:-1] for k in range(3)])
    xt = np.cumsum(dxt, axis=0) * dt
    xt -= xt[999:-1000].mean(0) - xn[999:-1000].mean(0)
    x, dx = xt[999:-1000], dxt[999:-1000]
    fig, ax = plt.subplots(3, 1, sharex=True)
    for k in range(3):
        ax[k].plot(dxt[:, k], "k", lw=0.8)
        ax[k].plot(dxclean[:, k], "r", lw=0.8)
        ax[k].set_xlim(5000, 7500)
else:
    dx = dxclean + 1.0 * rng.standard_normal(x.shape)
Xi = stls(pool_data(x, polyorder), dx, 0.025)
print_model(Xi, ["x", "y", "z"], polyorder)

model = sparse_galerkin(Xi, polyorder)
for T, tol in ((20, 1e-12), (250, 1e-6)):
    t = np.linspace(0, T, int(T / (0.001 if T == 20 else 0.01)) + 1)
    xA = simulate(lorenz, x0, t, rtol=tol, atol=tol)
    xB = simulate(model, x0, t, rtol=tol, atol=tol)
    fig = plt.figure(figsize=(10, 4))
    for k, xx in enumerate((xA, xB)):
        ax = fig.add_subplot(1, 2, k + 1, projection="3d")
        color_line3(ax, *xx.T, np.r_[0, np.diff(t)] + t / T)
        ax.view_init(16, 27)
        ax.set(xlabel="x", ylabel="y", zlabel="z", title=("true", "identified")[k])
    if T == 20:
        fig, ax = plt.subplots(1, 2, figsize=(10, 3))
        for k in range(2):
            ax[k].plot(t, xA[:, k], "k", t, xB[:, k], "r--", lw=1.5)
            ax[k].set(xlabel="Time", ylabel="xy"[k])
        print(f"t <= 20: trajectories stay within {np.abs(xA - xB)[t <= 5].max():.3g} up to t = 5")
plt.show()
