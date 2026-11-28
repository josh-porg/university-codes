"""Small fully connected sigmoid network trained by backpropagation (Math 796 ``Basic_NN``).

Port of Higham & Higham's ``activate``/``netbp``/``netbpfull``/``nlsrun``
("Deep learning: an introduction for applied mathematicians", 2019) and the
course's ``netbpfull_modified_*`` variants (vector packing for line searches
and black-box optimisers).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# The 10-point two-class data set of the paper: columns are points
HIGHAM_X = np.array([[0.1, 0.3, 0.1, 0.6, 0.4, 0.6, 0.5, 0.9, 0.4, 0.7],
                     [0.1, 0.4, 0.5, 0.9, 0.2, 0.3, 0.6, 0.2, 0.4, 0.6]])
HIGHAM_Y = np.vstack([np.r_[np.ones(5), np.zeros(5)], np.r_[np.zeros(5), np.ones(5)]])


def sigmoid(z):
    return 1 / (1 + np.exp(-z))


def activate(x, W, b):
    """``sigmoid(W x + b)`` (``activate``)."""
    return sigmoid(W @ x + b)


@dataclass
class SigmoidNetwork:
    """Network with layer widths ``sizes`` (e.g. ``(2, 2, 3, 2)``); ``x`` has one point per column."""

    sizes: tuple = (2, 2, 3, 2)
    weights: list = field(default_factory=list)
    biases: list = field(default_factory=list)

    @classmethod
    def random(cls, sizes=(2, 2, 3, 2), scale=0.5, rng=None):
        rng = np.random.default_rng(rng)
        W = [scale * rng.standard_normal((m, n)) for n, m in zip(sizes[:-1], sizes[1:])]
        b = [scale * rng.standard_normal((m, 1)) for m in sizes[1:]]
        return cls(tuple(sizes), W, b)

    def forward(self, x):
        """List of activations ``[x, a2, ..., aL]``."""
        a = [np.asarray(x, float)]
        for W, b in zip(self.weights, self.biases):
            a.append(activate(a[-1], W, b))
        return a

    def __call__(self, x):
        return self.forward(x)[-1]

    def residuals(self, x, y):
        """``|y_i - a_L(x_i)|`` per point (``neterr``); the cost is the sum of squares."""
        return np.linalg.norm(y - self(x), axis=0)

    def cost(self, x, y):
        return float(np.sum(self.residuals(x, y) ** 2))

    def gradients(self, x, y):
        """Backpropagated gradients of ``0.5 |a_L - y|^2`` summed over the columns of ``x``."""
        a = self.forward(x)
        delta = a[-1] * (1 - a[-1]) * (a[-1] - y)
        dW, db = [], []
        for layer in range(len(self.weights) - 1, -1, -1):
            dW.insert(0, delta @ a[layer].T)
            db.insert(0, delta.sum(axis=1, keepdims=True))
            if layer:
                delta = a[layer] * (1 - a[layer]) * (self.weights[layer].T @ delta)
        return dW, db

    def sgd(self, x, y, eta=0.05, iterations=100_000, rng=None, record_every=1000):
        """Stochastic gradient descent, one random point per step (``netbp``). Returns the cost history."""
        rng = np.random.default_rng(rng)
        history = []
        for k in range(iterations):
            i = rng.integers(x.shape[1])
            dW, db = self.gradients(x[:, i:i + 1], y[:, i:i + 1])
            for j in range(len(self.weights)):
                self.weights[j] -= eta * dW[j]
                self.biases[j] -= eta * db[j]
            if k % record_every == 0:
                history.append(self.cost(x, y))
        return np.array(history)

    # Vector packing (``vector_cost`` in the modified scripts). MATLAB order:
    # all weights (column-major) then all biases.
    def pack(self):
        return np.concatenate([W.ravel(order="F") for W in self.weights] + [b.ravel() for b in self.biases])

    def unpack(self, p):
        p = np.asarray(p, float)
        i = 0
        W, b = [], []
        for n, m in zip(self.sizes[:-1], self.sizes[1:]):
            W.append(p[i:i + m * n].reshape((m, n), order="F"))
            i += m * n
        for m in self.sizes[1:]:
            b.append(p[i:i + m].reshape(m, 1))
            i += m
        return SigmoidNetwork(self.sizes, W, b)

    def classify_grid(self, n=200):
        """Grid over the unit square and a mask where output 1 exceeds output 2 (the shaded decision regions)."""
        X, Y = np.meshgrid(np.linspace(0, 1, n + 1), np.linspace(0, 1, n + 1))
        out = self(np.vstack([X.ravel(), Y.ravel()]))
        return X, Y, (out[0] > out[1]).reshape(X.shape)
