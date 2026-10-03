"""Metrics of a structured 2D quadrilateral mesh."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class StructuredMesh2D:
    """Cell volumes and face geometry from node coordinates ``nodes[i, j] = (x, y)``.

    ``i_faces[i, j]`` is the face between cells ``i-1`` and ``i`` (shape
    ``(ni+1, nj, 3)``) and ``j_faces`` likewise along ``j``. Each holds
    ``[length, n_x, n_y]`` with the normal pointing towards increasing index.
    """

    nodes: np.ndarray

    def __post_init__(self):
        X = np.asarray(self.nodes, dtype=float)
        self.nodes = X
        # Cell areas from the cross product of the diagonals
        d1 = X[1:, 1:] - X[:-1, :-1]
        d2 = X[:-1, 1:] - X[1:, :-1]
        self.volume = 0.5 * np.abs(d1[..., 0] * d2[..., 1] - d1[..., 1] * d2[..., 0])
        self.centroid = 0.25 * (X[:-1, :-1] + X[1:, :-1] + X[1:, 1:] + X[:-1, 1:])
        self.i_faces = self._faces(X[:, 1:] - X[:, :-1], sign=1)  # edges along j
        self.j_faces = self._faces(X[1:, :] - X[:-1, :], sign=-1)  # edges along i

    @staticmethod
    def _faces(edge, sign):
        length = np.linalg.norm(edge, axis=-1)
        # rotate the edge by -90 deg (i faces) or +90 deg (j faces) to get the normal
        normal = sign * np.stack([edge[..., 1], -edge[..., 0]], axis=-1) / length[..., None]
        return np.concatenate([length[..., None], normal], axis=-1)

    @property
    def shape(self) -> tuple[int, int]:
        return self.volume.shape
