"""Finite-volume building blocks for the compressible Euler equations (AE 746).

State arrays keep the conserved variables on the *last* axis:
1D ``Q = [rho, rho u, E]``, 2D ``Q = [rho, rho u, rho v, E]``.

Face geometry for a structured 2D mesh is stored like the MATLAB
``AreaVI`` / ``AreaVJ``: ``[..., 0]`` face length, ``[..., 1:3]`` unit normal.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class FlowState:
    rho: np.ndarray
    velocity: np.ndarray  # (..., ndim)
    p: np.ndarray
    E: np.ndarray  # total energy per unit volume
    c: np.ndarray  # speed of sound

    @property
    def mach(self) -> np.ndarray:
        return np.linalg.norm(self.velocity, axis=-1) / self.c


def primitives(Q: np.ndarray, gamma: float = 1.4) -> FlowState:
    """Density, velocity, pressure, energy and sound speed from conserved ``Q``."""
    Q = np.asarray(Q, dtype=float)
    rho = Q[..., 0]
    mom = Q[..., 1:-1]
    E = Q[..., -1]
    vel = mom / rho[..., None]
    p = (gamma - 1) * (E - 0.5 * np.sum(mom * vel, axis=-1))
    c = np.sqrt(gamma * p / rho)
    return FlowState(rho, vel, p, E, c)


def conserved(rho, velocity, p, gamma: float = 1.4) -> np.ndarray:
    """Conserved state from density, velocity (..., ndim) and pressure."""
    rho = np.asarray(rho, dtype=float)
    velocity = np.asarray(velocity, dtype=float)
    E = p / (gamma - 1) + 0.5 * rho * np.sum(velocity**2, axis=-1)
    return np.concatenate([rho[..., None], rho[..., None] * velocity, np.asarray(E)[..., None]], axis=-1)


def normal_flux(Q, normal, gamma: float = 1.4) -> np.ndarray:
    """Physical flux of ``Q`` through a face with unit ``normal`` (..., ndim)."""
    s = primitives(Q, gamma)
    Vn = np.sum(s.velocity * normal, axis=-1)
    F = np.empty_like(np.asarray(Q, dtype=float))
    F[..., 0] = s.rho * Vn
    F[..., 1:-1] = s.rho[..., None] * s.velocity * Vn[..., None] + s.p[..., None] * normal
    F[..., -1] = Vn * (s.E + s.p)
    return F


def rusanov_flux(Q_L, Q_R, normal, gamma: float = 1.4) -> np.ndarray:
    """Rusanov (local Lax-Friedrichs) numerical flux across a face."""
    sL, sR = primitives(Q_L, gamma), primitives(Q_R, gamma)
    Vn_avg = 0.5 * (np.sum(sL.velocity * normal, -1) + np.sum(sR.velocity * normal, -1))
    c_avg = 0.5 * (sL.c + sR.c)
    dissipation = 0.5 * (np.abs(Vn_avg) + c_avg)[..., None] * (np.asarray(Q_R) - np.asarray(Q_L))
    return 0.5 * (normal_flux(Q_L, normal, gamma) + normal_flux(Q_R, normal, gamma)) - dissipation


def minmod_slope(Q, axis: int) -> np.ndarray:
    """Minmod-limited slope of interior cells along ``axis`` (length n-2)."""
    Q = np.moveaxis(np.asarray(Q, dtype=float), axis, 0)
    back = Q[1:-1] - Q[:-2]
    fwd = Q[2:] - Q[1:-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(back != 0, fwd / back, 0.0)
    slope = np.clip(r, 0.0, 1.0) * back
    return np.moveaxis(slope, 0, axis)


def reconstruct(Q, axis: int, order: int = 2):
    """Left and right states at the interior faces along ``axis``.

    For ``n`` cells returns arrays for the ``n-1`` interior faces. Cells
    next to the boundary fall back to first order, as in the MATLAB.
    """
    Q = np.moveaxis(np.asarray(Q, dtype=float), axis, 0)
    left = Q[:-1].copy()  # state on the low side of each face
    right = Q[1:].copy()  # state on the high side
    if order == 2:
        slope = np.moveaxis(minmod_slope(np.moveaxis(Q, 0, axis), axis), axis, 0)
        left[1:] = Q[1:-1] + 0.5 * slope
        right[:-1] = Q[1:-1] - 0.5 * slope
    elif order != 1:
        raise ValueError(f"order {order} not supported")
    return np.moveaxis(left, 0, axis), np.moveaxis(right, 0, axis)


# Ghost-state boundary conditions: each takes the interior state next to
# the boundary and returns the state on the far side of the face.


def slip_wall(Q_interior, normal, gamma: float = 1.4):
    """Inviscid wall: mirror the normal velocity component."""
    s = primitives(Q_interior, gamma)
    Vn = np.sum(s.velocity * normal, axis=-1, keepdims=True)
    return conserved(s.rho, s.velocity - 2 * Vn * normal, s.p, gamma)


def symmetry(Q_interior, gamma: float = 1.4, axis_component: int = 1):
    """Symmetry plane: flip the velocity component normal to the plane (default ``v``)."""
    Q_B = np.array(Q_interior, dtype=float, copy=True)
    Q_B[..., 1 + axis_component] *= -1
    return Q_B


def supersonic_exit(Q_interior, *_):
    """Supersonic outflow: extrapolate the interior state."""
    return np.array(Q_interior, dtype=float, copy=True)


def supersonic_inlet(Q_interior, Q_in):
    """Supersonic inflow: impose the free-stream state."""
    return np.broadcast_to(np.asarray(Q_in, dtype=float), np.shape(Q_interior)).copy()
