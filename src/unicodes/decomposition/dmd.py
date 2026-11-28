"""Exact dynamic mode decomposition (Tu et al. 2014; Kutz et al. 2016)."""

from __future__ import annotations

import numpy as np

from .linalg import svht_rank
from .snapshots import data_to_snapshots, snapshots_to_data


class DMD:
    """Exact DMD of a sequence of snapshots.

    Parameters
    ----------
    data : snapshot matrix ``(n_space, n_time)``, or gridded data with time
        as the last axis, e.g. ``(nx, ny, nt)`` or ``(nx, ny, nz, nt)``.
    rank : maximum number of modes to keep; ``None`` keeps all that pass
        the truncation.
    dt : time step between snapshots.
    truncation : how to discard negligible singular values before ``rank``
        is applied. A float drops singular values below that fraction of
        their sum (the MATLAB default was 0.005); ``"svht"`` uses the
        Gavish-Donoho optimal hard threshold; ``None`` keeps all.
    background_threshold : if given, modes with ``|omega|`` below this are
        treated as the slowly varying background and the rest as foreground
        (see :attr:`background` / :attr:`foreground`).

    Attributes
    ----------
    modes : DMD modes Phi, shape ``(n_space, r)``.
    eigenvalues : discrete-time eigenvalues lambda.
    omega : continuous-time eigenvalues ``log(lambda) / dt``.
    amplitudes : mode amplitudes b fitted to the first snapshot.
    singular_values : all singular values of the first snapshot matrix.
    rank : number of modes actually kept.
    """

    def __init__(
        self,
        data,
        rank: int | None = None,
        dt: float = 1.0,
        truncation: float | str | None = 0.005,
        background_threshold: float | None = None,
    ):
        data = np.asarray(data)
        if data.ndim == 2:
            self.spatial_shape: tuple[int, ...] | None = None
            self.snapshots = data
        else:
            self.spatial_shape = data.shape[:-1]
            self.snapshots = data_to_snapshots(data)

        self.dt = dt
        self.n_time = self.snapshots.shape[1]

        X, Xp = self.snapshots[:, :-1], self.snapshots[:, 1:]
        U, S, Vh = np.linalg.svd(X, full_matrices=False)
        self.singular_values = S

        if truncation is None:
            keep = len(S)
        elif truncation == "svht":
            keep = svht_rank(S, X.shape)
        else:
            keep = int(np.sum(S > truncation * S.sum()))
        r = keep if rank is None else min(rank, keep)
        if r < 1:
            raise ValueError("No singular values survive the truncation")
        self.rank = r

        U, S, V = U[:, :r], S[:r], Vh[:r].conj().T
        A_tilde = U.conj().T @ Xp @ V / S
        self.eigenvalues, W = np.linalg.eig(A_tilde)
        self._A_tilde, self._W, self._S, self._V = A_tilde, W, S, V
        self.modes = Xp @ V / S @ W  # exact DMD modes
        self.omega = np.log(self.eigenvalues.astype(complex)) / dt
        self.amplitudes = np.linalg.pinv(self.modes) @ self.snapshots[:, 0]

        self.background_threshold = background_threshold

    @property
    def frequencies(self) -> np.ndarray:
        """Mode frequencies in cycles per unit time, ``Im(omega) / 2 pi``."""
        return self.omega.imag / (2 * np.pi)

    def growth_adjusted_power(self) -> np.ndarray:
        """Mode power from amplitudes fitted to every snapshot (Gluhareff ``DMD``).

        ``b_t = (A_tilde W)^-1 Sigma V^*`` gives each mode's amplitude at every
        snapshot; the power is ``2 |b_t|`` summed in the 2-norm over time.
        Unlike ``2 |b|`` from the first snapshot alone, this is not biased by
        strongly growing or decaying modes.
        """
        b_t = np.linalg.solve(self._A_tilde @ self._W, self._S[:, None] * self._V.conj().T)
        return np.linalg.norm(2 * np.abs(b_t), axis=1)

    @property
    def time(self) -> np.ndarray:
        """Snapshot times starting at 0."""
        return np.arange(self.n_time) * self.dt

    def time_dynamics(self, t=None) -> np.ndarray:
        """``b * exp(omega t)`` for each mode, shape ``(rank, len(t))``."""
        t = self.time if t is None else np.asarray(t)
        return self.amplitudes[:, None] * np.exp(np.outer(self.omega, t))

    def reconstruct(self, t=None, modes=None) -> np.ndarray:
        """Reconstructed snapshots ``Phi exp(Omega t) b`` (complex).

        ``t`` defaults to the snapshot times, so this can also forecast.
        ``modes`` selects a subset of mode indices.
        """
        idx = slice(None) if modes is None else modes
        return self.modes[:, idx] @ self.time_dynamics(t)[idx]

    @property
    def reconstruction(self) -> np.ndarray:
        """Real part of the reconstruction, in the shape of the input data."""
        return self._to_data(self.reconstruct().real)

    def _split(self):
        if self.background_threshold is None:
            raise ValueError("Pass background_threshold to separate background and foreground")
        bg = np.abs(self.omega) < self.background_threshold
        return np.flatnonzero(bg), np.flatnonzero(~bg)

    @property
    def background(self) -> np.ndarray:
        """Reconstruction from the slow (|omega| < threshold) modes only."""
        return self._to_data(self._fit_subset(self._split()[0]).real)

    @property
    def foreground(self) -> np.ndarray:
        """Reconstruction from all other modes."""
        return self._to_data(self._fit_subset(self._split()[1]).real)

    def _fit_subset(self, idx):
        # Uses the amplitudes from the full fit. The MATLAB version refitted b
        # per subset, which leaks foreground content into the background
        # because DMD modes are not orthogonal.
        return self.reconstruct(modes=idx)

    def get_modes(self) -> np.ndarray:
        """Modes in the spatial shape of the input data, mode index last."""
        return self._to_data(self.modes)

    def _to_data(self, snapshots):
        if self.spatial_shape is None:
            return snapshots
        return snapshots_to_data(snapshots, self.spatial_shape)


def exact_dmd(X, Xp, rank=None, threshold=None):
    """Exact DMD of snapshot pairs ``Xp ~ A X`` that need not be a time sequence.

    Used for the MATH 526 diagnosis data, where each column pair is a
    patient's first and next hospital visit. Keeps singular values above
    ``threshold`` and at most ``rank`` of them. Returns
    ``(eigenvalues, modes, singular_values)``. (The MATLAB scaled the modes
    by the eigenvalues instead of using the eigenvectors.)
    """
    U, S, Vh = np.linalg.svd(np.asarray(X, float), full_matrices=False)
    r = len(S) if threshold is None else int(np.sum(S > threshold))
    r = r if rank is None else min(r, rank)
    U, Sr, V = U[:, :r], S[:r], Vh[:r].conj().T
    lam, W = np.linalg.eig(U.conj().T @ Xp @ V / Sr)
    return lam, Xp @ V / Sr @ W, S
