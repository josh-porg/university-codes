"""Dynamic mode decomposition with control (Proctor, Brunton & Kutz 2016; ``DMDc``)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class DMDcResult:
    A_tilde: np.ndarray  # reduced state matrix
    B_tilde: np.ndarray  # reduced input matrix
    U_hat: np.ndarray  # output-space basis: x ~ U_hat x_tilde
    modes: np.ndarray
    eigenvalues: np.ndarray
    omega: np.ndarray

    @property
    def A(self):
        return self.U_hat @ self.A_tilde @ self.U_hat.conj().T

    @property
    def B(self):
        return self.U_hat @ self.B_tilde


def _rank(S, tol):
    return int(np.sum(S > tol)) if tol is not None else len(S)


def dmdc(X, U, dt=1.0, tol=1e-5):
    """DMDc of states ``X`` (n, m) driven by inputs ``U`` (q, m).

    Singular values below ``tol`` are truncated in both SVDs (``None`` keeps
    all). The MATLAB took the input block of the left singular vectors as
    the single row ``n+q``; all ``q`` input rows are used here.
    """
    X, U = np.atleast_2d(X), np.atleast_2d(U)
    X0, X1, U0 = X[:, :-1], X[:, 1:], U[:, :-1]
    n = X.shape[0]
    Ut, St, Vth = np.linalg.svd(np.vstack([X0, U0]), full_matrices=False)
    p = _rank(St, tol)
    Ut, St, Vt = Ut[:, :p], St[:p], Vth[:p].conj().T
    Uh, Sh, _ = np.linalg.svd(X1, full_matrices=False)
    r = _rank(Sh, tol)
    Uh = Uh[:, :r]
    U1, U2 = Ut[:n], Ut[n:]
    core = X1 @ Vt / St
    A_t = Uh.conj().T @ core @ U1.conj().T @ Uh
    B_t = Uh.conj().T @ core @ U2.conj().T
    lam, W = np.linalg.eig(A_t)
    modes = core @ U1.conj().T @ Uh @ W
    return DMDcResult(A_t, B_t, Uh, modes, lam, np.log(lam.astype(complex)) / dt)
