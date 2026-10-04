"""Math 796: train Higham & Higham's 2-2-3-2 sigmoid network (``netbp``, ``netbpfull``, ``nlsrun``,
``netbpfull_modified_step``, ``netbpfull_modified_other_optimizers``).

Methods: ``sgd`` (stochastic gradient, eta = 0.05), ``lsq`` (nonlinear least
squares, ``lsqnonlin`` -> ``scipy.optimize.least_squares``), ``ga`` (genetic
algorithm on the packed parameter vector). ``--modified-labels`` uses the
relabelled targets of the modified scripts.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

from unicodes.neural import HIGHAM_X, HIGHAM_Y, SigmoidNetwork
from unicodes.optimize import genetic_minimize

p = argparse.ArgumentParser()
p.add_argument("--method", choices=["sgd", "lsq", "ga"], default="sgd")
p.add_argument("--iterations", type=int, default=1_000_000)
p.add_argument("--modified-labels", action="store_true")
args = p.parse_args()

x, y = HIGHAM_X, HIGHAM_Y
if args.modified_labels:
    first = np.array([0, 1, 1, 1, 1], bool)
    row = np.r_[first, ~first].astype(float)
    y = np.vstack([row, 1 - row])

net = SigmoidNetwork.random(rng=5000)  # MATLAB rng(5000); the random streams differ
if args.method == "sgd":
    history = net.sgd(x, y, eta=0.05, iterations=args.iterations, rng=0, record_every=10_000)
    fig, ax = plt.subplots()
    ax.semilogy(np.arange(history.size) * 10_000, history, "b-", lw=2)
    ax.set(xlabel="Iteration number", ylabel="Value of cost function")
elif args.method == "lsq":
    p0 = net.pack()
    r = least_squares(lambda q: net.unpack(q).residuals(x, y), p0)
    net = net.unpack(r.x)
else:
    p0 = net.pack()
    r = genetic_minimize(lambda q: net.unpack(q).cost(x, y), p0 - 5, p0 + 5, population_size=300, generations=500,
                         sigma=0.5, n_elites=10, rng=0)
    net = net.unpack(r.x)
print(f"{args.method}: final cost {net.cost(x, y):.4g}")

X, Y, mask = net.classify_grid(300)
fig, ax = plt.subplots()
ax.contourf(X, Y, mask.astype(float), [-0.5, 0.5, 1.5], colors=["white", "0.8"])
a = y[0] == 1
ax.plot(x[0, a], x[1, a], "ro", ms=12, mew=4, mfc="none")
ax.plot(x[0, ~a], x[1, ~a], "bx", ms=12, mew=4)
ax.set(xlim=(0, 1), ylim=(0, 1), xticks=[0, 1], yticks=[0, 1], title=f"Decision regions ({args.method})")
plt.show()
