"""Linear-algebra helpers for reduced-order modelling."""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy import integrate, optimize


def rsvd(X, rank: int, oversampling: int = 10, power_iterations: int = 1, rng=None):
    """Randomized truncated SVD of a (tall, skinny) matrix ``X``.

    Approximates the rank-``rank`` SVD. Oversampling by 5-10 is
    recommended; power iterations sharpen the singular-value spectrum at
    extra cost. Returns ``(U, S, Vh)`` like ``np.linalg.svd`` with
    ``full_matrices=False``, truncated to ``rank``.
    """
    rng = np.random.default_rng(rng)
    X = np.asarray(X)
    P = rng.standard_normal((X.shape[1], rank + oversampling))
    Z = X @ P
    for _ in range(power_iterations):
        Z = X @ (X.conj().T @ Z)
    Q, _ = np.linalg.qr(Z)
    Uy, S, Vh = np.linalg.svd(Q.conj().T @ X, full_matrices=False)
    return (Q @ Uy)[:, :rank], S[:rank], Vh[:rank]


def _lambda_star(beta):
    return np.sqrt(2 * (beta + 1) + 8 * beta / ((beta + 1) + np.sqrt(beta**2 + 14 * beta + 1)))


@lru_cache(maxsize=None)
def _marchenko_pastur_median(beta: float) -> float:
    lo, hi = (1 - np.sqrt(beta)) ** 2, (1 + np.sqrt(beta)) ** 2

    def density(t):
        return np.sqrt(max((hi - t) * (t - lo), 0.0)) / (2 * np.pi * beta * t)

    def cdf_minus_half(x):
        return integrate.quad(density, lo, x)[0] - 0.5

    return optimize.brentq(cdf_minus_half, lo, hi)


def optimal_svht_coef(beta, sigma_known: bool = False):
    """Optimal singular-value hard-threshold coefficient (Gavish & Donoho 2014).

    ``beta`` is the aspect ratio ``m/n <= 1`` of the matrix. With known
    noise level ``sigma`` keep singular values above
    ``coef * sqrt(n) * sigma``; with unknown noise keep those above
    ``coef * median(singular_values)``.

    Reference: M. Gavish and D. L. Donoho, "The Optimal Hard Threshold for
    Singular Values is 4/sqrt(3)", IEEE Trans. Inf. Theory 60(8), 2014.
    """
    beta_arr = np.atleast_1d(np.asarray(beta, dtype=float))
    if np.any((beta_arr <= 0) | (beta_arr > 1)):
        raise ValueError("beta must be in (0, 1]")
    coef = _lambda_star(beta_arr)
    if not sigma_known:
        coef = coef / np.sqrt([_marchenko_pastur_median(float(b)) for b in beta_arr])
    return coef if np.ndim(beta) else float(coef[0])


def svht_rank(singular_values, shape) -> int:
    """Number of singular values kept by the optimal hard threshold (unknown noise)."""
    s = np.asarray(singular_values)
    m, n = sorted(shape)
    return int(np.sum(s > optimal_svht_coef(m / n) * np.median(s)))


def cosamp(Phi, u, K, tol=1e-10, max_iterations=100):
    """Compressive sampling matching pursuit (Needell & Tropp 2009; book ``cosamp``).

    Finds a ``K``-sparse ``s`` with ``Phi s ~ u``. Works with complex ``Phi``.
    """
    Phi = np.asarray(Phi)
    u = np.asarray(u).ravel()
    s = np.zeros(Phi.shape[1], dtype=np.result_type(Phi, u))
    v = u.copy()
    support = np.array([], dtype=int)
    norm_u = np.linalg.norm(u)
    for _ in range(max_iterations):
        if np.linalg.norm(v) / norm_u <= tol:
            break
        y = np.abs(Phi.conj().T @ v)
        omega = np.flatnonzero((y >= np.sort(y)[::-1][min(2 * K, len(y)) - 1]) & (y > 1e-12))
        support = np.union1d(omega, support).astype(int)
        b = np.linalg.pinv(Phi[:, support]) @ u
        keep = np.abs(b) >= np.sort(np.abs(b))[::-1][min(K, len(b)) - 1]
        keep &= np.abs(b) > 1e-12
        support, b = support[keep], b[keep]
        s = np.zeros_like(s)
        s[support] = b
        v = u - Phi[:, support] @ b
    return s


def fast_ica(X, n_components, max_iterations=500, tol=1e-6, rng=None):
    """Independent components by symmetric FastICA with the tanh contrast (Hyvarinen 1999).

    ``X`` has one signal per row (signals, samples). Returns
    ``(sources, mixing)`` with ``sources`` of shape ``(n_components, samples)``
    and ``X - mean ~ mixing @ sources`` (the book's ``fastica``).
    """
    rng = np.random.default_rng(rng)
    X = np.asarray(X, float)
    Xc = X - X.mean(axis=1, keepdims=True)
    d, E = np.linalg.eigh(np.cov(Xc))
    idx = np.argsort(d)[::-1][:n_components]
    d, E = d[idx], E[:, idx]
    whiten = (E / np.sqrt(d)).T
    Z = whiten @ Xc
    W = np.linalg.qr(rng.standard_normal((n_components, n_components)))[0]
    for _ in range(max_iterations):
        G = np.tanh(W @ Z)
        W_new = G @ Z.T / Z.shape[1] - np.diag((1 - G**2).mean(axis=1)) @ W
        s, u = np.linalg.eigh(W_new @ W_new.T)
        W_new = u @ np.diag(1 / np.sqrt(s)) @ u.T @ W_new
        converged = np.max(np.abs(np.abs(np.diag(W_new @ W.T)) - 1)) < tol
        W = W_new
        if converged:
            break
    sources = W @ Z
    mixing = np.linalg.pinv(W @ whiten)
    return sources, mixing
