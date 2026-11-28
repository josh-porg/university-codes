"""Proper orthogonal decomposition (direct method, via the SVD)."""

from __future__ import annotations

import numpy as np

from .snapshots import data_to_snapshots, snapshots_to_data


class POD:
    """POD of snapshot data with the temporal mean removed.

    Parameters
    ----------
    data : snapshot matrix ``(n_space, n_time)`` or gridded data with time
        as the last axis.

    Attributes
    ----------
    mean : temporal mean, shape ``(n_space, 1)``.
    spatial_modes : orthonormal modes Phi, shape ``(n_space, k)``, ordered by energy.
    time_coefficients : a(t) for each mode, shape ``(n_time, k)``.
    eigenvalues : covariance eigenvalues ``sigma**2 / n_time`` (mode energies).

    The MATLAB version formed the ``n_space x n_space`` covariance matrix
    and called ``eig``; the SVD gives the same modes without building it.
    """

    def __init__(self, data):
        data = np.asarray(data)
        if data.ndim == 2:
            self.spatial_shape: tuple[int, ...] | None = None
            snapshots = data
        else:
            self.spatial_shape = data.shape[:-1]
            snapshots = data_to_snapshots(data)

        self.n_time = snapshots.shape[1]
        self.mean = snapshots.mean(axis=1, keepdims=True)
        fluct = snapshots - self.mean

        U, S, Vh = np.linalg.svd(fluct, full_matrices=False)
        self.spatial_modes = U
        self.time_coefficients = Vh.conj().T * S
        self.eigenvalues = S**2 / self.n_time

    @property
    def energy_fraction(self) -> np.ndarray:
        """Fraction of fluctuating energy captured by each mode."""
        return self.eigenvalues / self.eigenvalues.sum()

    def reconstruct(self, modes=None, include_mean: bool = True) -> np.ndarray:
        """Rebuild the data from a subset of modes (int ``k`` = first k, or indices)."""
        if modes is None:
            idx = slice(None)
        elif np.isscalar(modes):
            idx = slice(0, int(modes))
        else:
            idx = modes
        out = self.spatial_modes[:, idx] @ self.time_coefficients[:, idx].T
        if include_mean:
            out = out + self.mean
        return self._to_data(out)

    def get_modes(self) -> np.ndarray:
        """Spatial modes in the shape of the input data, mode index last."""
        return self._to_data(self.spatial_modes)

    def _to_data(self, snapshots):
        if self.spatial_shape is None:
            return snapshots
        return snapshots_to_data(snapshots, self.spatial_shape)
