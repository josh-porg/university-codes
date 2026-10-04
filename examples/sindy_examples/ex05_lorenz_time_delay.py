"""SINDy example 5 (``EX05_LorenzTimeDelay``): SINDy in time-delay coordinates of the Lorenz x signal.

The x component (dt = 0.001, t in (0, 100]) is shift-stacked 10 times; the
first three right singular vectors (eigen-time-delay coordinates u, v, w)
are differentiated by fourth-order central differences and a cubic
library with normalised columns is regressed with a different threshold
for each coordinate (0.01, 0.2, 2), as in the paper. The MATLAB allocated
ten Hankel rows but filled only the first three (the others stay zero and
do not change the SVD), and the thresholds are tuned for that, so the
three-delay embedding is kept (``--stacks`` changes it).
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from common import central_difference, pool_data, print_model, simulate, stls
from unicodes.decomposition import lorenz, time_delay_stack

p = argparse.ArgumentParser()
p.add_argument("--stacks", type=int, default=3)
a = p.parse_args()

dt = 0.001
tspan = np.arange(1, 100001) * dt
x = simulate(lorenz, [-8, 8, 27], tspan, rtol=1e-12, atol=1e-12)
stackmax, r = 10, 3
H = time_delay_stack(x[:len(x) - stackmax - 1 + a.stacks, 0], a.stacks)  # x(k:end-stackmax-1+k), k = 1..stacks
U, S, Vh = np.linalg.svd(H, full_matrices=False)
V = Vh[:r].T
dV = central_difference(V, dt)
X = V[2:-2]
Theta = pool_data(X, 3)
norms = np.linalg.norm(Theta, axis=0)
Xi = np.column_stack([stls(Theta / norms, dV[:, [k]], lam)[:, 0] for k, lam in enumerate((0.01, 0.2, 2))]) / norms[:, None]
print_model(Xi, ["u", "v", "w"], 3)
t = tspan[:len(X)]
fig, ax = plt.subplots(3, 1, sharex=True, figsize=(8, 6))
for k, name in enumerate("uvw"):
    ax[k].plot(t, Theta @ Xi[:, k], "k", label=r"$\Theta(V)\Xi$")
    ax[k].plot(t, dV[:, k], "r--", label=r"$\dot V$")
    ax[k].set(xlim=(40, 60), ylabel=name)
ax[0].legend()
ax[2].set_xlabel("Time")
plt.show()
