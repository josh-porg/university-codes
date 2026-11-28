"""Compressible Navier-Stokes terms for the structured 2-D solver: viscous fluxes, no-slip walls and
time-accurate integration (direct numerical simulation, no turbulence model).

Non-dimensional ideal gas ``p = rho R T``. The viscosity follows a power
law ``mu = mu_ref (T / T_ref)^omega`` (or Sutherland's law), conductivity
``k = mu c_p / Pr``. :func:`Viscosity.from_reynolds` sets ``mu_ref`` from a
free-stream Reynolds number ``rho U L / mu``.

DNS here means every scale is resolved by the grid and nothing is modelled;
the grid must resolve the boundary layers (``ruled_mesh(..., stretch=...)``
clusters points at the wall). In 2-D and at Reynolds numbers a few
thousand cells can resolve, the flows are laminar or transitional; DNS of
turbulence needs 3-D grids far beyond this solver.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import euler
from .mesh import StructuredMesh2D


@dataclass
class Viscosity:
    mu_ref: float
    T_ref: float = 1.0
    omega: float = 0.76  # power-law exponent; ignored when ``sutherland`` is given
    sutherland: float | None = None  # Sutherland temperature in the same units as T
    prandtl: float = 0.72
    R: float = 1.0

    @classmethod
    def from_reynolds(cls, reynolds, rho, speed, length, T_ref, **kw):
        return cls(mu_ref=rho * speed * length / reynolds, T_ref=T_ref, **kw)

    def mu(self, T):
        T = np.asarray(T, float)
        if self.sutherland is not None:
            S = self.sutherland
            return self.mu_ref * (T / self.T_ref) ** 1.5 * (self.T_ref + S) / (T + S)
        return self.mu_ref * (T / self.T_ref) ** self.omega

    def cp(self, gamma):
        return gamma * self.R / (gamma - 1)


def _face_geometry(mesh: StructuredMesh2D):
    X = mesh.nodes
    return 0.5 * (X[:, :-1] + X[:, 1:]), 0.5 * (X[:-1, :] + X[1:, :])  # i-face and j-face centres


def _primitives(Q, gamma, R):
    s = euler.primitives(Q, gamma)
    return np.stack([s.velocity[..., 0], s.velocity[..., 1], s.p / (s.rho * R)], axis=-1)  # u, v, T


def viscous_residual(mesh: StructuredMesh2D, Q, ghosts, visc: Viscosity, gamma=1.4):
    """Viscous contribution to dQ/dt.

    ``ghosts`` maps each boundary (``"i_min"`` ...) to its ghost conserved
    state (one layer, shape of the boundary cells), as built by the solver's
    boundary conditions. Gradients use Green-Gauss at cell centres, averaged
    to faces with the usual correction along the line joining the centroids.
    """
    phi = _primitives(Q, gamma, visc.R)
    xc = mesh.centroid
    fc_i, fc_j = _face_geometry(mesh)
    ni, nj = mesh.volume.shape

    g = {k: _primitives(v, gamma, visc.R) for k, v in ghosts.items()}
    # face values for Green-Gauss
    phi_i = np.concatenate([0.5 * (g["i_min"] + phi[:1]), 0.5 * (phi[1:] + phi[:-1]), 0.5 * (phi[-1:] + g["i_max"])], axis=0)
    phi_j = np.concatenate([0.5 * (g["j_min"] + phi[:, :1]), 0.5 * (phi[:, 1:] + phi[:, :-1]), 0.5 * (phi[:, -1:] + g["j_max"])], axis=1)
    Si, Ni = mesh.i_faces[..., 0], mesh.i_faces[..., 1:]
    Sj, Nj = mesh.j_faces[..., 0], mesh.j_faces[..., 1:]
    flux_i = phi_i[..., :, None] * (Ni * Si[..., None])[..., None, :]  # (ni+1, nj, 3 vars, 2 dims)
    flux_j = phi_j[..., :, None] * (Nj * Sj[..., None])[..., None, :]
    grad = (np.diff(flux_i, axis=0) + np.diff(flux_j, axis=1)) / mesh.volume[..., None, None]

    def face_gradient(phiL, phiR, gL, gR, xL, xR):
        avg = 0.5 * (gL + gR)
        dx = xR - xL
        d = np.linalg.norm(dx, axis=-1, keepdims=True)
        e = dx / d
        corr = (phiR - phiL) / d - np.einsum("...vk,...k->...v", avg, e)
        return avg + corr[..., None] * e[..., None, :]

    def mirror(xcell, xface, n):
        return xcell + 2 * np.sum((xface - xcell) * n, -1, keepdims=True) * n

    # i faces
    xL_i = np.concatenate([mirror(xc[:1], fc_i[:1], Ni[:1]), xc], axis=0)
    xR_i = np.concatenate([xc, mirror(xc[-1:], fc_i[-1:], Ni[-1:])], axis=0)
    pL_i = np.concatenate([g["i_min"], phi], axis=0)
    pR_i = np.concatenate([phi, g["i_max"]], axis=0)
    gL_i = np.concatenate([grad[:1], grad], axis=0)
    gR_i = np.concatenate([grad, grad[-1:]], axis=0)
    gf_i = face_gradient(pL_i, pR_i, gL_i, gR_i, xL_i, xR_i)
    # j faces
    xL_j = np.concatenate([mirror(xc[:, :1], fc_j[:, :1], Nj[:, :1]), xc], axis=1)
    xR_j = np.concatenate([xc, mirror(xc[:, -1:], fc_j[:, -1:], Nj[:, -1:])], axis=1)
    pL_j = np.concatenate([g["j_min"], phi], axis=1)
    pR_j = np.concatenate([phi, g["j_max"]], axis=1)
    gL_j = np.concatenate([grad[:, :1], grad], axis=1)
    gR_j = np.concatenate([grad, grad[:, -1:]], axis=1)
    gf_j = face_gradient(pL_j, pR_j, gL_j, gR_j, xL_j, xR_j)

    k_factor = visc.cp(gamma) / visc.prandtl

    def flux(pL, pR, gf, N, S):
        u, v, T = (0.5 * (pL[..., k] + pR[..., k]) for k in range(3))
        mu = visc.mu(T)
        ux, uy, vx, vy, Tx, Ty = gf[..., 0, 0], gf[..., 0, 1], gf[..., 1, 0], gf[..., 1, 1], gf[..., 2, 0], gf[..., 2, 1]
        div = ux + vy
        txx = mu * (2 * ux - 2 / 3 * div)
        tyy = mu * (2 * vy - 2 / 3 * div)
        txy = mu * (uy + vx)
        kx, ky = mu * k_factor * Tx, mu * k_factor * Ty
        nx, ny = N[..., 0], N[..., 1]
        F = np.stack([np.zeros_like(u), txx * nx + txy * ny, txy * nx + tyy * ny,
                      (u * txx + v * txy + kx) * nx + (u * txy + v * tyy + ky) * ny], axis=-1)
        return F * S[..., None]

    Fi = flux(pL_i, pR_i, gf_i, Ni, Si)
    Fj = flux(pL_j, pR_j, gf_j, Nj, Sj)
    assert Fi.shape[:2] == (ni + 1, nj) and Fj.shape[:2] == (ni, nj + 1)
    return (np.diff(Fi, axis=0) + np.diff(Fj, axis=1)) / mesh.volume[..., None]


def viscous_time_step_denominator(mesh: StructuredMesh2D, Q, visc: Viscosity, gamma=1.4):
    """Viscous part of the time-step denominator, ``max(4/3, gamma/Pr) (mu/rho) sum(S^2) / V``."""
    s = euler.primitives(Q, gamma)
    nu = visc.mu(s.p / (s.rho * visc.R)) / s.rho
    S2 = (mesh.i_faces[:-1, :, 0] ** 2 + mesh.i_faces[1:, :, 0] ** 2 + mesh.j_faces[:, :-1, 0] ** 2
          + mesh.j_faces[:, 1:, 0] ** 2)
    return 2 * max(4 / 3, gamma / visc.prandtl) * nu * S2 / mesh.volume


@dataclass
class UnsteadyResult:
    Q: np.ndarray
    time: float
    steps: int
    snapshots: list = field(default_factory=list)  # (t, Q) pairs
    monitor: list = field(default_factory=list)  # (t, value) from the monitor callback
