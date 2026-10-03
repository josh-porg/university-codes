"""Steady 2D Euler solver on a structured mesh (AE 746 project 4).

Rusanov fluxes, optional MUSCL/minmod reconstruction, SSP-RK2 pseudo-time
stepping with local (or global) time steps.

    bcs = {"i_min": "symmetry", "i_max": "exit", "j_min": "wall", "j_max": ("inlet", Q_in)}
    result = solve_steady(mesh, Q0, bcs, order=2, cfl=0.5)
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import euler
from .mesh import StructuredMesh2D


@dataclass
class SteadyResult:
    Q: np.ndarray
    residual_history: list[float] = field(default_factory=list)
    converged: bool = False
    iterations: int = 0


def _ghost(kind, Q_edge, faces, gamma):
    if isinstance(kind, tuple):
        kind, Q_in = kind
    if kind == "wall":
        return euler.slip_wall(Q_edge, faces[..., 1:], gamma)
    if kind == "symmetry":
        # mirror about the face: equivalent to a slip wall
        return euler.slip_wall(Q_edge, faces[..., 1:], gamma)
    if kind == "exit":
        return euler.supersonic_exit(Q_edge)
    if kind == "inlet":
        return euler.supersonic_inlet(Q_edge, Q_in)
    raise ValueError(f"Unknown boundary condition {kind!r}")


def _face_states(Q, axis, faces, bc_min, bc_max, order, gamma):
    left, right = euler.reconstruct(Q, axis, order)
    first = np.take(Q, [0], axis=axis)
    last = np.take(Q, [-1], axis=axis)
    lo_faces = np.take(faces, [0], axis=axis)
    hi_faces = np.take(faces, [-1], axis=axis)
    ghost_lo = _ghost(bc_min, first, lo_faces, gamma)
    ghost_hi = _ghost(bc_max, last, hi_faces, gamma)
    L = np.concatenate([ghost_lo, left, last], axis=axis)
    R = np.concatenate([first, right, ghost_hi], axis=axis)
    return L, R


def residual(mesh: StructuredMesh2D, Q, bcs, order=2, gamma=1.4, return_faces=False):
    """dQ/dt for every cell (flux balance divided by cell volume)."""
    out = []
    total = 0.0
    for axis, faces, lo, hi in ((0, mesh.i_faces, "i_min", "i_max"), (1, mesh.j_faces, "j_min", "j_max")):
        L, R = _face_states(Q, axis, faces, bcs[lo], bcs[hi], order, gamma)
        F = euler.rusanov_flux(L, R, faces[..., 1:], gamma) * faces[..., :1]
        net = np.diff(F, axis=axis)
        total = total + net
        out.append((L, R, faces))
    res = -total / mesh.volume[..., None]
    return (res, out) if return_faces else res


def local_time_step(mesh, face_data, cfl, gamma=1.4, global_step=False):
    denom = 0.0
    for axis, (L, R, faces) in enumerate(face_data):
        sL, sR = euler.primitives(L, gamma), euler.primitives(R, gamma)
        n = faces[..., 1:]
        Vn = 0.5 * (np.sum(sL.velocity * n, -1) + np.sum(sR.velocity * n, -1))
        wave = (np.abs(Vn) + 0.5 * (sL.c + sR.c)) * faces[..., 0]
        lo = np.take(wave, range(wave.shape[axis] - 1), axis=axis)
        hi = np.take(wave, range(1, wave.shape[axis]), axis=axis)
        denom = denom + lo + hi
    dt = 2 * cfl * mesh.volume / denom
    return np.full_like(dt, dt.min()) if global_step else dt


def solve_steady(
    mesh: StructuredMesh2D,
    Q0,
    bcs: dict,
    order: int = 2,
    cfl: float = 0.5,
    gamma: float = 1.4,
    tol: float = 1e-6,
    max_iterations: int = 20000,
    global_step: bool = False,
    callback=None,
) -> SteadyResult:
    """March to steady state with SSP-RK2 until the density residual drops by ``tol``.

    ``bcs`` maps ``"i_min"``, ``"i_max"``, ``"j_min"``, ``"j_max"`` to one of
    ``"wall"``, ``"symmetry"``, ``"exit"`` or ``("inlet", Q_in)``.
    """
    Q = np.array(Q0, dtype=float, copy=True)
    result = SteadyResult(Q)
    res0 = None
    for it in range(1, max_iterations + 1):
        res, faces = residual(mesh, Q, bcs, order, gamma, return_faces=True)
        dt = local_time_step(mesh, faces, cfl, gamma, global_step)[..., None]
        Q_star = Q + dt * res
        Q = 0.5 * (Q + Q_star + dt * residual(mesh, Q_star, bcs, order, gamma))

        norm = float(np.linalg.norm(res[..., 0]))
        res0 = res0 or norm or 1.0
        result.residual_history.append(norm)
        if callback is not None:
            callback(it, Q, norm)
        if norm <= res0 * tol:
            result.converged = True
            break
    result.Q, result.iterations = Q, it
    return result
