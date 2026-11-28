"""Flux-reconstruction (DG-recovering) operators on the reference element [-1, 1] (AE 846 ``FR_coefs``)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.polynomial import legendre as L


def lagrange_basis(nodes, x):
    """Values ``l_i(x)`` of the Lagrange basis on ``nodes``; shape ``(len(x), len(nodes))``."""
    nodes, x = np.asarray(nodes, float), np.atleast_1d(np.asarray(x, float))
    out = np.ones((x.size, nodes.size))
    for i, xi in enumerate(nodes):
        for j, xj in enumerate(nodes):
            if i != j:
                out[:, i] *= (x - xj) / (xi - xj)
    return out


def differentiation_matrix(nodes):
    """``D[i, j] = l_j'(x_i)``."""
    nodes = np.asarray(nodes, float)
    n = nodes.size
    D = np.zeros((n, n))
    for j in range(n):
        others = np.delete(nodes, j)
        c = np.atleast_1d(np.poly(others)) / np.prod(nodes[j] - others)
        D[:, j] = np.polyval(np.polyder(c), nodes)
    return D


@dataclass
class FROperators:
    nodes: np.ndarray  # Gauss-Legendre solution points
    D: np.ndarray  # differentiation matrix
    interp: np.ndarray  # rows: values at x = -1 and x = 1
    dg_left: np.ndarray  # derivative of the left Radau correction at the nodes
    dg_right: np.ndarray


def fr_operators(k):
    """Operators for degree-``k`` solution polynomials with the DG (right/left Radau) correction.

    ``g_L = (-1)^(k+1)/2 (P_(k+1) - P_k)``, ``g_R(x) = g_L(-x)``. The MATLAB
    reused ``k`` as a symbolic variable and always built the k = 2
    correction; any ``k`` works here.
    """
    nodes, _ = L.leggauss(k + 1)
    gL = (-1) ** (k + 1) / 2 * (L.Legendre.basis(k + 1) - L.Legendre.basis(k))
    dgl = gL.deriv()(nodes)
    return FROperators(nodes, differentiation_matrix(nodes), lagrange_basis(nodes, [-1.0, 1.0]), dgl, -dgl[::-1])
