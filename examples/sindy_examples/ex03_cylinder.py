"""SINDy example 3 (``EX03_Cylinder``): a mean-field model of the cylinder wake from POD coefficients.

Needs the paper's ``PODcoefficients.mat`` and ``PODcoefficients_run1.mat``
(``alpha``, ``alphaS``; not in this repository): the first two POD
coefficients and the shift mode of two DNS runs at dt = 0.02. Derivatives
by fourth-order central differences, fifth-order library, threshold 2e-5.
The identified model has constant terms because the fixed point is not at
the origin in these coordinates; it is integrated from the start of each
run::

    python ex03_cylinder.py PODcoefficients.mat PODcoefficients_run1.mat
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from common import central_difference, color_line3, pool_data, print_model, simulate, sparse_galerkin, stls
from unicodes.io import load_mat

p = argparse.ArgumentParser()
p.add_argument("run0")
p.add_argument("run1")
a = p.parse_args()

dt, r, polyorder = 0.02, 2, 5


def load(path, n):
    d = load_mat(path)
    alpha, alphaS = np.asarray(d["alpha"]), np.asarray(d["alphaS"]).reshape(len(d["alphaS"]), -1)
    x = np.column_stack([alpha[:n, :r], alphaS[:n, 0]])
    return x[2:-2], central_difference(x, dt)


x0_, dx0 = load(a.run0, 5000)
x1_, dx1 = load(a.run1, 3000)
x, dx = np.vstack([x0_, x1_]), np.vstack([dx0, dx1])
Xi = stls(pool_data(x, polyorder), dx, 2e-5)
print_model(Xi, ["x", "y", "z"], polyorder)
model = sparse_galerkin(Xi, polyorder)
for start, data, T in ((0, x0_, 100), (len(x0_), x1_, 95)):
    t = np.arange(0, T + dt / 2, dt)
    xD = simulate(model, x[start], t, rtol=1e-8, atol=1e-8)
    fig = plt.figure(figsize=(10, 4))
    for k, xx in enumerate((data, xD)):
        ax = fig.add_subplot(1, 2, k + 1, projection="3d")
        color_line3(ax, xx[:, 0], xx[:, 1], xx[:, -1], dt * np.arange(len(xx)))
        ax.view_init(16, 27)
        ax.set(xlim=(-200, 200), ylim=(-200, 200), zlim=(-160, 20), title=("DNS", "identified")[k])
plt.show()
