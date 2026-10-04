"""Quasi-1D Euler equations in a converging-diverging nozzle (AE 746 project 2).

Ports ``area``, ``MachArea``, ``Initialize``, ``ExactSolu`` and the finite-volume
solver of ``Project_2_Main_V_*``/``Project_2_MUSCL``: cell states on the
``x`` points, MUSCL reconstruction (none/minmod/van Leer/Barth limiters),
Rusanov fluxes, the source ``-(A'/A) [rho u, rho u^2, u (E + p)]``,
SSP-RK2 pseudo-time marching with global or local steps, and either fixed
(exact-solution) or characteristic boundary states. Non-dimensional with
``rho = 1``, ``p = 1/gamma`` at the reference station. ``Q`` is ``(n, 3)``.

Fixes from the MATLAB: the Rusanov dissipation used ``|u|`` (or each
eigenvalue separately) and was halved twice; the standard spectral radius
``|u| + c`` is used.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq

from . import euler

GAMMA = 1.4


def nozzle_area(x):
    """Course nozzle area distribution, throat (A = 0.338486) at x = 0 (``area``)."""
    x = np.asarray(x, dtype=float)
    g = np.exp(-np.log(2.0) * (x / 0.6) ** 2)
    return np.where(x >= 0, 0.536572 - 0.198086 * g, 1.0 - 0.661514 * g)


def nozzle_area_slope(x):
    x = np.asarray(x, dtype=float)
    dg = np.exp(-np.log(2.0) * (x / 0.6) ** 2) * (-2 * np.log(2.0) * x / 0.36)
    return np.where(x >= 0, -0.198086 * dg, -0.661514 * dg)


def _area_ratio(M, gamma=GAMMA):
    return 1 / M * (2 / (gamma + 1) * (1 + (gamma - 1) / 2 * M**2)) ** ((gamma + 1) / (2 * (gamma - 1)))


def _mach(A_ratio, supersonic, gamma=GAMMA):
    if A_ratio <= 1 + 1e-12:
        return 1.0
    f = lambda M: _area_ratio(M, gamma) - A_ratio  # noqa: E731
    return brentq(f, 1 + 1e-12, 50) if supersonic else brentq(f, 1e-6, 1 - 1e-12)


def _state(M, p0, rho0, gamma=GAMMA):
    t = 1 + 0.5 * (gamma - 1) * M**2
    p, rho = p0 / t ** (gamma / (gamma - 1)), rho0 / t ** (1 / (gamma - 1))
    u = M * np.sqrt(gamma * p / rho)
    return np.array([rho, rho * u, p / (gamma - 1) + 0.5 * rho * u**2])


@dataclass
class NozzleFlow:
    """Exact quasi-1D solution (``Initialize`` + ``ExactSolu``).

    ``case=1``: subsonic throughout with exit Mach 0.4 at x = 5.
    ``case=2``: inlet Mach 0.2006533, exit pressure 0.6071752 at x = 4, giving a
    normal shock in the diverging section (located by bisection as in the MATLAB).
    """

    case: int = 1
    gamma: float = GAMMA
    x_shock: float = np.inf
    p0: float = field(init=False)
    rho0: float = field(init=False)
    a_star: float = field(init=False)
    p02: float = field(init=False, default=np.nan)
    a_star2: float = field(init=False, default=np.nan)

    def __post_init__(self):
        g = self.gamma
        if self.case == 1:
            M_e = 0.4
            t = 1 + 0.5 * (g - 1) * M_e**2
            self.p0, self.rho0 = (1 / g) * t ** (g / (g - 1)), t ** (1 / (g - 1))
            self.a_star = float(nozzle_area(5.0)) / _area_ratio(M_e, g)
        elif self.case == 2:
            M_i, p_e = 0.2006533, 0.6071752
            t = 1 + 0.5 * (g - 1) * M_i**2
            self.p0, self.rho0 = (1 / g) * t ** (g / (g - 1)), t ** (1 / (g - 1))
            self.a_star = float(nozzle_area(0.0))
            self.x_shock = brentq(lambda xs: self._exit_pressure(xs) - p_e, 1e-3, 3.999)
            self._exit_pressure(self.x_shock)
        else:
            raise ValueError("case must be 1 or 2")

    def _exit_pressure(self, xs):
        g = self.gamma
        A = float(nozzle_area(xs))
        Ms = _mach(A / self.a_star, True, g)
        p1 = self.p0 / (1 + 0.5 * (g - 1) * Ms**2) ** (g / (g - 1))
        p2 = p1 * (1 + 2 * g / (g + 1) * (Ms**2 - 1))
        M2 = np.sqrt((1 + 0.5 * (g - 1) * Ms**2) / (g * Ms**2 - 0.5 * (g - 1)))
        self.a_star2 = A / _area_ratio(M2, g)
        self.p02 = p2 * (1 + 0.5 * (g - 1) * M2**2) ** (g / (g - 1))
        Me = _mach(float(nozzle_area(4.0)) / self.a_star2, False, g)
        return self.p02 / (1 + 0.5 * (g - 1) * Me**2) ** (g / (g - 1))

    def mach(self, x):
        x = float(x)
        if x > self.x_shock:
            return _mach(float(nozzle_area(x)) / self.a_star2, False, self.gamma)
        supersonic = self.case == 2 and x > 0
        return _mach(float(nozzle_area(x)) / self.a_star, supersonic, self.gamma)

    def __call__(self, x):
        """Conserved state(s) at ``x``."""
        xs = np.atleast_1d(np.asarray(x, dtype=float))
        out = []
        for xi in xs:
            M = self.mach(xi)
            if xi > self.x_shock:
                T0 = self.p0 / self.rho0
                rho0_2 = self.p02 / T0  # stagnation temperature is unchanged across the shock
                out.append(_state(M, self.p02, rho0_2, self.gamma))
            else:
                out.append(_state(M, self.p0, self.rho0, self.gamma))
        out = np.array(out)
        return out[0] if np.ndim(x) == 0 else out


def characteristic_state(Q_L, Q_R, gamma=GAMMA):
    """Boundary state from the entropy and Riemann invariants of the left and right states (``findBoundaryQ``)."""
    sL, sR = euler.primitives(Q_L, gamma), euler.primitives(Q_R, gamma)
    uL, uR = sL.velocity[..., 0], sR.velocity[..., 0]
    s = sL.p / sL.rho**gamma
    Jp = uL + 2 * sL.c / (gamma - 1)
    Jm = uR - 2 * sR.c / (gamma - 1)
    u = 0.5 * (Jp + Jm)
    c = (gamma - 1) / 4 * (Jp - Jm)
    rho = (c**2 / (gamma * s)) ** (1 / (gamma - 1))
    p = rho * c**2 / gamma
    return euler.conserved(rho, np.asarray(u)[..., None], p, gamma)


def _limited_slopes(Q, Q_in, Q_out, dx, limiter):
    Qe = np.vstack([Q_in, Q, Q_out])
    dL = (Qe[1:-1] - Qe[:-2]) / dx
    dR = (Qe[2:] - Qe[1:-1]) / dx
    if limiter in (None, "none"):
        return 0.5 * (dL + dR)
    if limiter == "minmod":
        return np.where(dL * dR > 0, np.sign(dL) * np.minimum(abs(dL), abs(dR)), 0.0)
    if limiter == "vanleer":
        return (np.sign(dL) + np.sign(dR)) * dL * dR / (abs(dL) + abs(dR) + 1e-30)
    if limiter == "barth":
        slope = 0.5 * (dL + dR)
        Qmin = np.minimum(Qe[1:-1], np.minimum(Qe[:-2], Qe[2:]))
        Qmax = np.maximum(Qe[1:-1], np.maximum(Qe[:-2], Qe[2:]))
        phi = np.ones_like(Q)
        for sign in (-1, 1):
            dQ = sign * 0.5 * dx * slope
            with np.errstate(divide="ignore", invalid="ignore"):
                lim = np.where(dQ > 0, (Qmax - Q) / dQ, np.where(dQ < 0, (Qmin - Q) / dQ, 1.0))
            phi = np.minimum(phi, np.clip(lim, 0, 1))
        return phi * slope
    raise ValueError(f"unknown limiter {limiter!r}")


def residual(Q, x, Q_in, Q_out, order=2, limiter=None, characteristic=True, gamma=GAMMA):
    """``dQ/dt`` for the cell states ``Q`` at points ``x`` (uniform spacing)."""
    dx = x[1] - x[0]
    if characteristic:
        Q_in = characteristic_state(Q_in, Q[0], gamma)
        Q_out = characteristic_state(Q[-1], Q_out, gamma)
    S = _limited_slopes(Q, Q_in, Q_out, dx, limiter) if order == 2 else np.zeros_like(Q)
    left_face = Q - 0.5 * dx * S  # state on the left face of each cell
    right_face = Q + 0.5 * dx * S
    QL = np.vstack([Q_in, right_face])  # faces 0..n: state on the left side
    QR = np.vstack([left_face, Q_out])
    F = euler.rusanov_flux(QL, QR, np.ones((len(QL), 1)), gamma)
    s = euler.primitives(Q, gamma)
    u = s.velocity[:, 0]
    k = -nozzle_area_slope(x) / nozzle_area(x)
    G = np.column_stack([k * s.rho * u, k * s.rho * u**2, k * u * (s.E + s.p)])
    return (F[:-1] - F[1:]) / dx + G


@dataclass
class NozzleResult:
    x: np.ndarray
    Q: np.ndarray
    residual_history: list
    converged: bool


def solve_nozzle(flow: NozzleFlow, n_cells=100, x_range=(-4.0, 4.0), order=2, limiter=None, cfl=0.25,
                 global_step=True, characteristic=True, tol=1e-6, max_iterations=200_000) -> NozzleResult:
    """March to steady state from the exit state; stops when |R| falls by ``tol``."""
    x = np.linspace(*x_range, n_cells)
    dx = x[1] - x[0]
    Q_in, Q_out = flow(x_range[0]), flow(x_range[1])
    Q = np.tile(Q_out, (n_cells, 1))
    g = flow.gamma
    R0 = None
    history = []
    for _ in range(max_iterations):
        s = euler.primitives(Q, g)
        dt = cfl * dx / np.abs(s.velocity[:, 0] + s.c)
        if global_step:
            dt = np.full_like(dt, dt.min())
        R = residual(Q, x, Q_in, Q_out, order, limiter, characteristic, g)
        Q_star = Q + dt[:, None] * R
        R_star = residual(Q_star, x, Q_in, Q_out, order, limiter, characteristic, g)
        Q = 0.5 * (Q + Q_star + dt[:, None] * R_star)
        norm = np.linalg.norm(R_star)
        R0 = norm if R0 is None else R0
        history.append(np.linalg.norm(R_star[:, 0]))
        if norm <= tol * R0:
            return NozzleResult(x, Q, history, True)
    return NozzleResult(x, Q, history, False)
