"""DMD variants from Kutz, Brunton, Brunton & Proctor, *Dynamic Mode Decomposition* (SIAM 2016).

Ports of the book's companion code: forward-backward and total-least-squares
DMD (``DMD_eig``, ch. 8), optimal mode amplitudes (Jovanovic et al. 2014,
used by ``mrDMD``), multiresolution DMD (``mrDMD``, ``mrDMD_map``, ch. 5),
the eigensystem realization algorithm (``ERA``, ch. 7), time-delay stacking
(ch. 7, 12) and kernel DMD (``Algorithm_10_7through11``, ch. 10).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import sqrtm


def _atilde(Y, Yp, r=None):
    U, S, Vh = np.linalg.svd(Y, full_matrices=False)
    r = len(S) if r is None else min(r, len(S))
    U, S, V = U[:, :r], S[:r], Vh[:r].conj().T
    return U.conj().T @ Yp @ V / S, U, S, V


def dmd_eigenvalues(data, flavor="dmd", rank=None):
    """Discrete-time DMD eigenvalues of snapshots ``data`` (n, m) (``DMD_eig``).

    ``flavor``: ``"dmd"`` (exact DMD), ``"fbdmd"`` (forward-backward DMD,
    ``(A_f A_b^-1)^(1/2)``, which removes the noise bias of DMD; both
    operators are formed in the POD basis of the first snapshots, as in
    Dawson et al. 2016 — the book code took the backward one in the basis of
    the shifted snapshots, which biases even noise-free data) or
    ``"tlsdmd"`` (total-least-squares DMD). The book's ``DMD_eig`` sized the
    TLS blocks by the number of snapshots; they are sized by the state
    dimension (or ``rank``) here.
    """
    Y, Yp = data[:, :-1], data[:, 1:]
    flavor = flavor.lower()
    if flavor == "dmd":
        A = _atilde(Y, Yp, rank)[0]
    elif flavor == "fbdmd":
        U = _atilde(Y, Yp, rank)[1]
        Xt, Yt = U.conj().T @ Y, U.conj().T @ Yp  # both directions in the same POD basis
        Af = Yt @ np.linalg.pinv(Xt)
        Ab = Xt @ np.linalg.pinv(Yt)
        A = sqrtm(Af @ np.linalg.inv(Ab))
    elif flavor == "tlsdmd":
        r = data.shape[0] if rank is None else rank
        if rank is not None and rank < data.shape[0]:  # project onto the leading POD subspace first
            Ux = np.linalg.svd(Y, full_matrices=False)[0][:, :r]
            Y, Yp = Ux.conj().T @ Y, Ux.conj().T @ Yp
        U = np.linalg.svd(np.vstack([Y, Yp]), full_matrices=False)[0]
        A = U[r:2 * r, :r] @ np.linalg.inv(U[:r, :r])
    else:
        raise ValueError(f"unknown DMD flavor {flavor!r}")
    return np.linalg.eigvals(A)


def optimal_amplitudes(eigenvalues, W, S, V):
    """Amplitudes ``b`` that best fit all snapshots (Jovanovic, Schmid & Nichols 2014).

    ``W`` are the eigenvectors of ``A_tilde`` and ``S``, ``V`` the truncated
    singular values and right singular vectors of the first snapshot matrix.
    """
    m = V.shape[0]
    vand = eigenvalues[:, None] ** np.arange(m)[None, :]
    G = S[:, None] * V.conj().T
    P = (W.conj().T @ W) * np.conj(vand @ vand.conj().T)
    q = np.conj(np.diag(vand @ G.conj().T @ W))
    L = np.linalg.cholesky(P)
    return np.linalg.solve(L.conj().T, np.linalg.solve(L, q))


def time_delay_stack(X, stacks):
    """Shift-stack snapshots ``X`` (n, m) ``stacks`` times: shape ``(n stacks, m - stacks + 1)``."""
    X = np.atleast_2d(X)
    m = X.shape[1]
    return np.vstack([X[:, k:m - stacks + 1 + k] for k in range(stacks)])


@dataclass
class MrDMDLevel:
    level: int
    start: int  # first snapshot of this window
    stop: int
    T: float  # window length
    rho: float  # frequency cutoff (cycles per unit time)
    omega: np.ndarray  # slow-mode frequencies, cycles per unit time (complex)
    power: np.ndarray  # |b| of the slow modes
    modes: np.ndarray  # slow modes (2 n rows: the data and one delay)


@dataclass
class MrDMD:
    levels: int
    n_time: int
    windows: list = field(default_factory=list)

    def level(self, k):
        return [w for w in self.windows if w.level == k]

    def amplitude_map(self, resolution=None):
        """Time-frequency map of mean slow-mode power (``mrDMD_map``) and the level cutoffs.

        Returns ``(map, low_f)`` with ``map`` of shape ``(levels, resolution)``
        (default ``2^(levels-1)`` columns) and ``low_f[k]`` the lower frequency
        bound of level ``k``.
        """
        M = resolution or 2 ** (self.levels - 1)
        out = np.zeros((self.levels, M))
        low_f = np.zeros(self.levels + 1)
        for k in range(self.levels):
            for w in self.level(k):
                P = w.power[np.abs(w.omega.imag) >= low_f[k]]
                if P.size:
                    a = int(round(w.start / self.n_time * M))
                    b = int(round(w.stop / self.n_time * M))
                    out[k, a:b] = P.mean()
            low_f[k + 1] = self.level(k)[0].rho
        return out, low_f


def mrdmd(X, dt, rank, max_cycles=2, levels=6) -> MrDMD:
    """Multiresolution DMD (Kutz, Fu & Brunton 2016; ``mrDMD``).

    At each level the window is subsampled to about four times the Nyquist
    rate of ``rho = max_cycles / T``, augmented with one delay, decomposed
    with rank ``rank`` DMD and optimal amplitudes; modes slower than
    ``rho`` are kept for that level and the two halves of the window are
    processed at the next level.
    """
    X = np.asarray(X)
    out = MrDMD(levels, X.shape[1])

    def recurse(start, stop, level):
        Xw = X[:, start:stop]
        T = Xw.shape[1] * dt
        rho = max_cycles / T
        sub = max(int(np.ceil(1 / rho / 8 / np.pi / dt)), 1)
        Xs = Xw[:, ::sub]
        Xaug = np.vstack([Xs[:, :-1], Xs[:, 1:]])
        A, U, S, V = _atilde(Xaug[:, :-1], Xaug[:, 1:], rank)
        lam, W = np.linalg.eig(A)
        Phi = Xaug[:, 1:] @ V / S @ W
        b = optimal_amplitudes(lam, W, S, V)
        omega = np.log(lam.astype(complex)) / sub / dt / (2 * np.pi)
        slow = np.abs(omega) <= rho
        out.windows.append(MrDMDLevel(level, start, stop, T, rho, omega[slow], np.abs(b[slow]), Phi[:, slow]))
        if level + 1 < levels:
            mid = start + (stop - start) // 2
            recurse(start, mid, level + 1)
            recurse(mid, stop, level + 1)

    recurse(0, X.shape[1], 0)
    return out


def era(Y, m, n, r):
    """Eigensystem realization algorithm (Juang & Pappa 1985; ``ERA``).

    ``Y`` holds impulse-response Markov parameters with shape
    shape (outputs, inputs, steps); ``Y[:, :, 0]`` is the feedthrough D.
    ``m`` x ``n`` block Hankel matrices, reduced order ``r``. Returns
    ``(A, B, C, D, hankel_singular_values)``.
    """
    Y = np.asarray(Y, float)
    nout, nin, _ = Y.shape
    D = Y[:, :, 0]
    Yk = Y[:, :, 1:]
    if Yk.shape[2] < m + n:
        raise ValueError("need at least m + n Markov parameters after the first")
    H = np.empty((nout * m, nin * n))
    H2 = np.empty_like(H)
    for i in range(m):
        for j in range(n):
            H[i * nout:(i + 1) * nout, j * nin:(j + 1) * nin] = Yk[:, :, i + j]
            H2[i * nout:(i + 1) * nout, j * nin:(j + 1) * nin] = Yk[:, :, i + j + 1]
    U, S, Vh = np.linalg.svd(H, full_matrices=False)
    Ur, Sr, Vr = U[:, :r], S[:r], Vh[:r].T
    s_half = 1 / np.sqrt(Sr)
    A = (s_half[:, None] * Ur.T) @ H2 @ (Vr * s_half[None, :])
    B = (s_half[:, None] * Ur.T) @ H[:, :nin]
    C = H[:nout, :] @ (Vr * s_half[None, :])
    return A, B, C, D, S


def kernel_dmd(X1, X2, kernel, rank=None, tol=1e-12):
    """Kernel DMD eigenvalues (Williams, Rowley & Kevrekidis 2015).

    ``kernel(a, b)`` is evaluated between snapshot columns: the Gram matrix
    ``G_ij = k(x_i, x_j)`` on ``X1`` and ``A_ij = k(y_i, x_j)`` with ``Y = X2``.
    With ``G = Q Sigma^2 Q^*`` the Koopman matrix ``K = (Sigma^+ Q^*) A (Q Sigma^+)``
    has the eigenvalues of ``A G^+``, computed here with a pseudo-inverse
    truncated at ``tol`` (relative) or to ``rank`` singular values, which also
    covers the non-Hermitian kernels of the book. The book code built the
    Gram matrix from the X1-X2 pairs too and took the eigenvalues of
    ``Sigma^2 Q^T A Q Sigma^-2`` (which are those of ``Q^T A Q``).
    """
    m = X1.shape[1]
    G = np.array([[kernel(X1[:, i], X1[:, j]) for j in range(m)] for i in range(m)])
    A = np.array([[kernel(X2[:, i], X1[:, j]) for j in range(m)] for i in range(m)])
    U, S, Vh = np.linalg.svd(G)
    keep = S > tol * S[0]
    if rank is not None:
        keep &= np.arange(len(S)) < rank
    # nonzero eigenvalues of A G^+ = A V S^-1 U^*: those of the r x r matrix U^* A V S^-1
    return np.linalg.eigvals(U[:, keep].conj().T @ A @ Vh[keep].conj().T / S[keep])
