"""Finite-difference schemes for the periodic 1-D linear advection equation ``u_t + c u_x = 0``
(AE 746 project 1, ``Project_1``/``Project_1_part_2``).

Each step function takes the periodic solution vector and the CFL number
``sigma = c dt / dx`` (c > 0).
"""

from __future__ import annotations

import numpy as np


def upwind1(u, sigma):
    """First-order upwind, forward Euler in time."""
    return u - sigma * (u - np.roll(u, 1))


def central2_euler(u, sigma):
    """Second-order central space, forward Euler time (unconditionally unstable)."""
    return u - sigma / 2 * (np.roll(u, -1) - np.roll(u, 1))


def upwind2_residual(u, sigma):
    """``sigma/2 (3 u_i - 4 u_{i-1} + u_{i-2})``."""
    return sigma / 2 * (3 * u - 4 * np.roll(u, 1) + np.roll(u, 2))


def upwind2_ssprk2(u, sigma):
    """Second-order upwind space, SSP-RK2 (Heun) time.

    The first ``Project_1`` applied the predictor twice; ``Project_1_part_2``
    (used here) averages ``u^n`` and the corrected predictor.
    """
    u_star = u - upwind2_residual(u, sigma)
    return 0.5 * (u + u_star - upwind2_residual(u_star, sigma))


SCHEMES = {"upwind1": upwind1, "central2_euler": central2_euler, "upwind2_ssprk2": upwind2_ssprk2}


def advect(u0, sigma, steps, scheme="upwind2_ssprk2"):
    """March ``steps`` time steps; returns the final solution."""
    step = SCHEMES[scheme] if isinstance(scheme, str) else scheme
    u = np.asarray(u0, dtype=float).copy()
    for _ in range(steps):
        u = step(u, sigma)
    return u
