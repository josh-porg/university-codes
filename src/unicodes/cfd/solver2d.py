"""2D Euler / Navier-Stokes solver on a structured mesh (AE 746 project 4).

Rusanov fluxes, optional MUSCL/minmod reconstruction of the primitive
variables, SSP-RK2 pseudo-time stepping with local (or global) time steps
for steady flow, and SSP-RK3 with a global time step for time-accurate
simulation. Passing a :class:`~unicodes.cfd.navier_stokes.Viscosity` adds
the viscous and heat-conduction fluxes (laminar Navier-Stokes / DNS).

    bcs = {"i_min": "symmetry", "i_max": "exit", "j_min": "wall", "j_max": ("inlet", Q_in)}
    result = solve_steady(mesh, Q0, bcs, order=2, cfl=0.5)

Boundary kinds: ``"wall"`` (inviscid slip), ``"noslip"`` (adiabatic no-slip
wall) or ``("noslip", T_wall)`` (isothermal, needs ``R`` of the viscosity
model), ``"symmetry"``, ``"exit"`` (supersonic extrapolation) and
``("inlet", Q_in)``.
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


def _ghost(kind, Q_edge, faces, gamma, R=1.0):
    if isinstance(kind, tuple):
        kind, Q_in = kind
    else:
        Q_in = None
    if kind == "noslip":
        st = euler.primitives(Q_edge, gamma)
        if Q_in is None:  # adiabatic: same temperature, reversed velocity
            return euler.conserved(st.rho, -st.velocity, st.p, gamma)
        T_cell = st.p / (st.rho * R)
        T_ghost = np.maximum(2 * Q_in - T_cell, 0.05 * Q_in)  # wall temperature is the face average
        return euler.conserved(st.p / (R * T_ghost), -st.velocity, st.p, gamma)
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


def _face_states(Q, axis, faces, bc_min, bc_max, order, gamma, R_gas=1.0):
    left, right = euler.reconstruct(Q, axis, order, gamma)
    first = np.take(Q, [0], axis=axis)
    last = np.take(Q, [-1], axis=axis)
    lo_faces = np.take(faces, [0], axis=axis)
    hi_faces = np.take(faces, [-1], axis=axis)
    ghost_lo = _ghost(bc_min, first, lo_faces, gamma, R_gas)
    ghost_hi = _ghost(bc_max, last, hi_faces, gamma, R_gas)
    L = np.concatenate([ghost_lo, left, last], axis=axis)
    R = np.concatenate([first, right, ghost_hi], axis=axis)
    return L, R, ghost_lo, ghost_hi


def residual(mesh: StructuredMesh2D, Q, bcs, order=2, gamma=1.4, return_faces=False, viscosity=None):
    """dQ/dt for every cell (flux balance divided by cell volume)."""
    out = []
    total = 0.0
    ghosts = {}
    R_gas = viscosity.R if viscosity is not None else 1.0
    for axis, faces, lo, hi in ((0, mesh.i_faces, "i_min", "i_max"), (1, mesh.j_faces, "j_min", "j_max")):
        L, R, g_lo, g_hi = _face_states(Q, axis, faces, bcs[lo], bcs[hi], order, gamma, R_gas)
        ghosts[lo], ghosts[hi] = g_lo, g_hi
        F = euler.rusanov_flux(L, R, faces[..., 1:], gamma) * faces[..., :1]
        net = np.diff(F, axis=axis)
        total = total + net
        out.append((L, R, faces))
    res = -total / mesh.volume[..., None]
    if viscosity is not None:
        from .navier_stokes import viscous_residual

        res = res + viscous_residual(mesh, Q, ghosts, viscosity, gamma)
    return (res, out) if return_faces else res


def local_time_step(mesh, face_data, cfl, gamma=1.4, global_step=False, Q=None, viscosity=None):
    denom = 0.0
    if viscosity is not None:
        from .navier_stokes import viscous_time_step_denominator

        denom = 2 * viscous_time_step_denominator(mesh, Q, viscosity, gamma)
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
    viscosity=None,
) -> SteadyResult:
    """March to steady state with SSP-RK2 until the density residual drops by ``tol``.

    ``bcs`` maps ``"i_min"``, ``"i_max"``, ``"j_min"``, ``"j_max"`` to one of
    ``"wall"``, ``"symmetry"``, ``"exit"`` or ``("inlet", Q_in)``.
    """
    Q = np.array(Q0, dtype=float, copy=True)
    result = SteadyResult(Q)
    res0 = None
    for it in range(1, max_iterations + 1):
        res, faces = residual(mesh, Q, bcs, order, gamma, return_faces=True, viscosity=viscosity)
        dt = local_time_step(mesh, faces, cfl, gamma, global_step, Q, viscosity)[..., None]
        Q_star = Q + dt * res
        Q = 0.5 * (Q + Q_star + dt * residual(mesh, Q_star, bcs, order, gamma, viscosity=viscosity))

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


def solve_unsteady(mesh: StructuredMesh2D, Q0, bcs: dict, t_end: float, order: int = 2, cfl: float = 0.5,
                   gamma: float = 1.4, viscosity=None, snapshot_every: float | None = None, monitor=None,
                   max_steps: int = 10_000_000):
    """Time-accurate integration to ``t_end`` with SSP-RK3 and a global time step.

    With a ``viscosity`` this is a direct numerical simulation of the
    compressible Navier-Stokes equations (nothing modelled). ``monitor(t, Q)``
    may return a number to record each step (e.g. a drag or a probe value);
    ``snapshot_every`` stores ``(t, Q)`` copies.
    """
    from .navier_stokes import UnsteadyResult

    Q = np.array(Q0, dtype=float, copy=True)
    t, k = 0.0, 0
    out = UnsteadyResult(Q, 0.0, 0)
    next_snap = 0.0

    def L(Qs):
        return residual(mesh, Qs, bcs, order, gamma, viscosity=viscosity)

    while t < t_end - 1e-14 and k < max_steps:
        res, faces = residual(mesh, Q, bcs, order, gamma, return_faces=True, viscosity=viscosity)
        dt = float(local_time_step(mesh, faces, cfl, gamma, True, Q, viscosity).min())
        dt = min(dt, t_end - t)
        Q1 = Q + dt * res
        Q2 = 0.75 * Q + 0.25 * (Q1 + dt * L(Q1))
        Q = Q / 3 + 2 / 3 * (Q2 + dt * L(Q2))
        t += dt
        k += 1
        if not np.all(np.isfinite(Q)):
            raise FloatingPointError(f"solution blew up at t = {t:.4g} (step {k}); reduce cfl")
        if monitor is not None:
            val = monitor(t, Q)
            if val is not None:
                out.monitor.append((t, val))
        if snapshot_every is not None and t >= next_snap - 1e-12:
            out.snapshots.append((t, Q.copy()))
            next_snap += snapshot_every
    out.Q, out.time, out.steps = Q, t, k
    return out
