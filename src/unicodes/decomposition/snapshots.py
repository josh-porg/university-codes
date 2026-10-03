"""Converting gridded data to and from snapshot matrices."""

from __future__ import annotations

import numpy as np


def data_to_snapshots(data: np.ndarray) -> np.ndarray:
    """Flatten data with time as the last axis into a snapshot matrix.

    ``(nx, ny, ..., nt)`` -> ``(nx*ny*..., nt)``: each column is one
    snapshot. Undo with :func:`snapshots_to_data`.
    """
    data = np.asarray(data)
    return data.reshape(-1, data.shape[-1])


def snapshots_to_data(snapshots: np.ndarray, spatial_shape: tuple[int, ...]) -> np.ndarray:
    """Reshape a snapshot matrix ``(n_space, nt)`` back to ``(*spatial_shape, nt)``."""
    snapshots = np.asarray(snapshots)
    return snapshots.reshape(*spatial_shape, snapshots.shape[-1])


def remove_mean(snapshots: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Subtract the temporal mean from each row of a snapshot matrix.

    Returns ``(mean_free_snapshots, mean)`` where ``mean`` has shape
    ``(n_space, 1)`` so it broadcasts back on with ``mean_free + mean``.
    """
    snapshots = np.asarray(snapshots)
    mean = snapshots.mean(axis=1, keepdims=True)
    return snapshots - mean, mean
