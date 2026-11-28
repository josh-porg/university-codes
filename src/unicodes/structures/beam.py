"""Euler-Bernoulli beam finite elements and stiffness sensitivities (AE 725 final exam problem 4)."""

from __future__ import annotations

import numpy as np


def beam_element_stiffness(E, I, L):
    """4x4 stiffness of a two-node beam element, DOFs ``[w1, theta1, w2, theta2]``."""
    return E * I / L**3 * np.array([
        [12, 6 * L, -12, 6 * L],
        [6 * L, 4 * L**2, -6 * L, 2 * L**2],
        [-12, -6 * L, 12, -6 * L],
        [6 * L, 2 * L**2, -6 * L, 4 * L**2],
    ], dtype=float)  # fmt: skip


def assemble_beam(E, I, L):
    """Global stiffness of a chain of beam elements with properties ``E[i], I[i], L[i]``."""
    E, I, L = np.broadcast_arrays(np.atleast_1d(E), np.atleast_1d(I), np.atleast_1d(L))
    n = E.size
    K = np.zeros((2 * n + 2, 2 * n + 2))
    for i in range(n):
        K[2 * i:2 * i + 4, 2 * i:2 * i + 4] += beam_element_stiffness(E[i], I[i], L[i])
    return K


def displacement_sensitivity(K, dK, F, free, output, adjoint=False):
    """Derivative of ``output . u`` with respect to a design variable whose stiffness derivative is ``dK``.

    ``free`` indexes the unconstrained DOFs. Direct: ``-p^T K^-1 dK u``;
    adjoint: ``lambda = K^-T p``, ``-lambda^T dK u`` (identical result, one
    solve per output instead of per variable).
    """
    Kf = K[np.ix_(free, free)]
    dKf = dK[np.ix_(free, free)]
    u = np.linalg.solve(Kf, np.asarray(F, float)[free])
    p = np.asarray(output, float)
    if adjoint:
        lam = np.linalg.solve(Kf.T, p)
        return -lam @ dKf @ u
    return -p @ np.linalg.solve(Kf, dKf @ u)
