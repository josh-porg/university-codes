"""Small numerical helpers shared by several ports."""

from __future__ import annotations

import numpy as np


def jacobian(f, x, step=1e-3, relative=True):
    """Central-difference Jacobian of ``f: R^n -> R^m`` at ``x`` (``computeJacobian``, AE 722/725).

    ``relative=True`` perturbs each variable by ``step * |x_j|`` (falling
    back to ``step`` where x_j = 0); otherwise by ``step``. The MATLAB
    "fixed step" option actually used ``|x_j| + step``, which is not a fixed
    step; use ``relative=False`` for a true fixed step.
    """
    x = np.asarray(x, dtype=float)
    f0 = np.atleast_1d(f(x))
    J = np.zeros((f0.size, x.size))
    for j in range(x.size):
        h = step * abs(x[j]) if relative and x[j] != 0 else step
        e = np.zeros_like(x)
        e[j] = h
        J[:, j] = (np.atleast_1d(f(x + e)) - np.atleast_1d(f(x - e))) / (2 * h)
    return J


def gradient(f, x, step=1e-6):
    """Central-difference gradient of a scalar function."""
    return jacobian(lambda z: np.atleast_1d(f(z)), x, step, relative=False)[0]


def smoothed_noise(n, windows=(50, 5), rng=None):
    """White noise smoothed by successive moving averages (``random_smoothed_number_playground``)."""
    rng = np.random.default_rng(rng)
    x = rng.standard_normal(n)
    for w in windows:
        x = np.convolve(x, np.ones(w) / w)[:n]  # causal like MATLAB filter()
    return x
