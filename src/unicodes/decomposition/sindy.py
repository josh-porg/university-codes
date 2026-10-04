"""Sparse identification of nonlinear dynamics (SINDy; Brunton, Proctor & Kutz 2016).

Ports ``generateLibrary``, ``generateLibraryList``, ``STLS``, ``STRidge``,
``InterpretResults`` and ``SINDyGalerkin`` from the DMD/SINDy project.
States are rows: ``X`` has shape ``(n_samples, n_states)``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations_with_replacement

import numpy as np
from scipy.integrate import solve_ivp

POLY_ORDERS = {"const": 0, "linear": 1, "poly2": 2, "poly3": 3, "poly4": 4, "poly5": 5}


def library(X, terms=("const", "linear", "poly2"), n_harmonics=10):
    """Candidate-function matrix Theta(X) (``generateLibrary``).

    ``terms`` from ``"const"``, ``"linear"``, ``"poly2"`` ... ``"poly5"`` (all
    monomials of that degree) and ``"sin"`` (``sin(k x)``, ``cos(k x)`` for
    k = 1..n_harmonics).
    """
    X = np.atleast_2d(np.asarray(X, dtype=float))
    cols = []
    n = X.shape[1]
    for term in terms:
        if term in POLY_ORDERS:
            for idx in combinations_with_replacement(range(n), POLY_ORDERS[term]):
                cols.append(np.prod(X[:, list(idx)], axis=1) if idx else np.ones(len(X)))
        elif term == "sin":
            for k in range(1, n_harmonics + 1):
                cols.extend(np.sin(k * X).T)
                cols.extend(np.cos(k * X).T)
        else:
            raise ValueError(f"unknown library term {term!r}")
    return np.column_stack(cols)


def library_names(names, terms=("const", "linear", "poly2"), n_harmonics=10):
    """Labels for the columns of :func:`library` (``generateLibraryList``)."""
    out = []
    for term in terms:
        if term in POLY_ORDERS:
            for idx in combinations_with_replacement(range(len(names)), POLY_ORDERS[term]):
                out.append("".join(names[i] for i in idx) if idx else "1")
        elif term == "sin":
            for k in range(1, n_harmonics + 1):
                out.extend(f"sin({k}{v})" for v in names)
                out.extend(f"cos({k}{v})" for v in names)
    return out


def stls(Theta, dXdt, threshold, iterations=10):
    """Sequentially thresholded least squares (``STLS``)."""
    Xi = np.linalg.lstsq(Theta, dXdt, rcond=None)[0]
    for _ in range(iterations):
        small = np.abs(Xi) < threshold
        Xi[small] = 0
        for j in range(dXdt.shape[1]):
            big = ~small[:, j]
            if big.any():
                Xi[big, j] = np.linalg.lstsq(Theta[:, big], dXdt[:, j], rcond=None)[0]
    return Xi


def stridge(Theta, dXdt, threshold, ridge=0.0, iterations=10):
    """Sequentially thresholded ridge regression (``STRidge``); equals :func:`stls` for ``ridge=0``."""
    n_terms = Theta.shape[1]
    Xi = np.zeros((n_terms, dXdt.shape[1]))
    for j in range(dXdt.shape[1]):
        big = np.ones(n_terms, bool)
        xi = np.zeros(n_terms)
        for _ in range(iterations + 1):
            T = Theta[:, big]
            xi[:] = 0
            xi[big] = np.linalg.solve(T.T @ T + ridge * np.eye(big.sum()), T.T @ dXdt[:, j]) if big.any() else 0
            new_big = np.abs(xi) >= threshold
            if np.array_equal(new_big, big):
                break
            big = new_big
        xi[~big] = 0
        Xi[:, j] = xi
    return Xi


@dataclass
class SINDy:
    """Fit ``dX/dt = Theta(X) Xi`` and simulate the identified model."""

    terms: tuple = ("const", "linear", "poly2")
    threshold: float = 0.1
    ridge: float = 0.0
    n_harmonics: int = 10
    coefficients: np.ndarray | None = field(default=None, init=False)

    def fit(self, X, dXdt):
        Theta = library(X, self.terms, self.n_harmonics)
        self.coefficients = stridge(Theta, np.asarray(dXdt, float), self.threshold, self.ridge)
        return self

    def predict(self, X):
        """Model derivatives at the states ``X``."""
        return library(X, self.terms, self.n_harmonics) @ self.coefficients

    def rhs(self, t, x):
        """Right-hand side for ODE solvers (``SINDyGalerkin``)."""
        return self.predict(np.atleast_2d(x))[0]

    def simulate(self, x0, t, **kwargs):
        kwargs.setdefault("rtol", 1e-10)
        kwargs.setdefault("atol", 1e-10)
        sol = solve_ivp(self.rhs, (t[0], t[-1]), x0, t_eval=t, **kwargs)
        return sol.y.T

    def equations(self, names, precision=4):
        """Readable identified equations (``InterpretResults``)."""
        labels = library_names(names, self.terms, self.n_harmonics)
        out = []
        for j, name in enumerate(names):
            parts = [f"{c:+.{precision}g} {lab}" for c, lab in zip(self.coefficients[:, j], labels) if c != 0]
            out.append(f"d{name}/dt = " + (" ".join(parts) if parts else "0"))
        return out


def lorenz(t, x, sigma=10.0, rho=28.0, beta=8 / 3):
    """Lorenz system (``lorenz``)."""
    return np.array([sigma * (x[1] - x[0]), x[0] * (rho - x[2]) - x[1], x[0] * x[1] - beta * x[2]])


def van_der_pol(t, x, mu=1.0):
    """Van der Pol oscillator as a first-order system (``vanderpol``)."""
    return np.array([x[1], mu * (1 - x[0] ** 2) * x[1] - x[0]])
