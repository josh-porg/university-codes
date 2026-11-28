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


def tv_derivative(data, iterations, alpha, dx=None, ep=1e-6, scale="small", u0=None):
    """Total-variation regularised derivative of noisy data (``TVRegDiff``).

    Port of R. Chartrand's ``TVRegDiff`` ("Numerical differentiation of noisy,
    nonsmooth data", ISRN Applied Mathematics, 2011): lagged-diffusivity
    iterations with conjugate-gradient inner solves. ``scale="small"``
    returns ``len(data) + 1`` values on the cell edges (the SINDy scripts
    dropped the last one); ``"large"`` returns ``len(data)`` values and suits
    long signals (the original returned it per sample; it is divided by
    ``dx`` here so both variants are derivatives).
    """
    from scipy.sparse import diags
    from scipy.sparse.linalg import LinearOperator, cg

    data = np.asarray(data, dtype=float).ravel()
    n = data.size
    dx = 1 / n if dx is None else dx
    if scale == "small":
        D = diags([-np.ones(n), np.ones(n)], [0, 1], shape=(n, n + 1)) / dx

        def A(x):
            return (np.cumsum(x) - 0.5 * (x + x[0]))[1:] * dx

        def AT(w):
            return (w.sum() - np.r_[w.sum() / 2, np.cumsum(w) - w / 2]) * dx

        u = np.r_[0, np.diff(data), 0] if u0 is None else np.asarray(u0, float)
        ATb = AT(data[0] - data)
        for _ in range(iterations):
            Q = diags(1 / np.sqrt((D @ u) ** 2 + ep))
            Lm = dx * (D.T @ Q @ D)
            g = AT(A(u)) + ATb + alpha * (Lm @ u)
            op = LinearOperator((n + 1, n + 1), matvec=lambda v: alpha * (Lm @ v) + AT(A(v)))
            P = LinearOperator((n + 1, n + 1), matvec=lambda v, d=alpha * (Lm.diagonal() + 1): v / d)
            s, _ = cg(op, g, rtol=1e-4, maxiter=100, M=P)
            u = u - s
        return u
    if scale == "large":
        D = diags([-np.ones(n), np.ones(n - 1)], [0, 1], shape=(n, n)).tolil() / dx
        D[n - 1, n - 1] = 0
        D = D.tocsr()

        def A(v):
            return np.cumsum(v)

        def AT(w):
            return w.sum() - np.r_[0, np.cumsum(w[:-1])]

        d = data - data[0]
        u = np.r_[0, np.diff(d)] if u0 is None else np.asarray(u0, float)
        ATd = AT(d)
        c = np.cumsum(np.arange(n, 0, -1))[::-1]
        for _ in range(iterations):
            Q = diags(1 / np.sqrt((D @ u) ** 2 + ep))
            Lm = D.T @ Q @ D
            g = AT(A(u)) - ATd + alpha * (Lm @ u)
            op = LinearOperator((n, n), matvec=lambda v: alpha * (Lm @ v) + AT(A(v)))
            P = LinearOperator((n, n), matvec=lambda v, dd=alpha * Lm.diagonal() + c: v / dd)
            s, _ = cg(op, -g, rtol=1e-4, maxiter=100, M=P)
            u = u + s
        return u / dx  # Chartrand's "large" variant integrates without dx, so u is per sample
    raise ValueError("scale must be 'small' or 'large'")
