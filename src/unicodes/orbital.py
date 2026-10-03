"""Two-body orbital mechanics (AE 360).

SI units (m, m/s, s, rad). ``mu`` defaults to Earth's.

    r, v = propagate_kepler(r0, v0, dt=86400)          # analytic
    sol = propagate(r0, v0, t_span=(0, 86400))           # numerical (solve_ivp)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

MU_EARTH = 3.986004418e14  # m^3/s^2


def two_body_rhs(t, x, mu=MU_EARTH):
    """State derivative for ``x = [r, v]``."""
    r = x[:3]
    return np.concatenate([x[3:], -mu * r / np.linalg.norm(r) ** 3])


def specific_energy(r, v, mu=MU_EARTH):
    return 0.5 * np.sum(np.square(v), axis=-1) - mu / np.linalg.norm(r, axis=-1)


def specific_angular_momentum(r, v):
    return np.cross(r, v)


def propagate_euler(r0, v0, duration, dt, mu=MU_EARTH):
    """Explicit-Euler propagation (as in the homework); returns ``(t, r, v)`` histories.

    Kept for teaching comparisons: energy drifts with step size. Prefer :func:`propagate`.
    """
    n = int(round(duration / dt))
    x = np.concatenate([r0, v0]).astype(float)
    out = np.empty((n + 1, 6))
    out[0] = x
    for k in range(n):
        x = x + dt * two_body_rhs(0, x, mu)
        out[k + 1] = x
    t = np.arange(n + 1) * dt
    return t, out[:, :3], out[:, 3:]


def propagate(r0, v0, t_span, mu=MU_EARTH, rtol=1e-10, atol=1e-6, **kwargs):
    """Numerically propagate with ``scipy.integrate.solve_ivp`` (DOP853)."""
    return solve_ivp(
        two_body_rhs, t_span, np.concatenate([r0, v0]).astype(float), args=(mu,),
        method="DOP853", rtol=rtol, atol=atol, **kwargs,
    )


def solve_kepler(M, e, tol=1e-12, max_iter=50):
    """Eccentric anomaly ``E`` from mean anomaly ``M`` (Newton's method)."""
    M = np.mod(M + np.pi, 2 * np.pi) - np.pi  # wrap to (-pi, pi]
    E = M - e if -np.pi < M < 0 else M + e
    for _ in range(max_iter):
        dE = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        E -= dE
        if abs(dE) < tol:
            return E
    raise RuntimeError("Kepler's equation did not converge")


@dataclass(frozen=True)
class OrbitalElements:
    a: float  # semi-major axis (m)
    e: float  # eccentricity
    i: float  # inclination (rad)
    raan: float  # right ascension of the ascending node (rad)
    argp: float  # argument of periapsis (rad)
    nu: float  # true anomaly (rad)


def rv_to_coe(r, v, mu=MU_EARTH) -> OrbitalElements:
    """Classical orbital elements from position and velocity (elliptic, non-equatorial orbits)."""
    r, v = np.asarray(r, float), np.asarray(v, float)
    rn, vn = np.linalg.norm(r), np.linalg.norm(v)
    h = np.cross(r, v)
    hn = np.linalg.norm(h)
    n = np.cross([0, 0, 1], h)
    nn = np.linalg.norm(n)
    e_vec = ((vn**2 - mu / rn) * r - np.dot(r, v) * v) / mu
    e = np.linalg.norm(e_vec)
    a = 1 / (2 / rn - vn**2 / mu)
    i = np.arccos(h[2] / hn)
    raan = np.arccos(np.clip(n[0] / nn, -1, 1))
    if n[1] < 0:
        raan = 2 * np.pi - raan
    argp = np.arccos(np.clip(np.dot(n, e_vec) / (nn * e), -1, 1))
    if e_vec[2] < 0:
        argp = 2 * np.pi - argp
    nu = np.arccos(np.clip(np.dot(e_vec, r) / (e * rn), -1, 1))
    if np.dot(r, v) < 0:
        nu = 2 * np.pi - nu
    return OrbitalElements(a, e, i, raan, argp, nu)


def coe_to_rv(coe: OrbitalElements, mu=MU_EARTH):
    """Position and velocity from classical orbital elements."""
    a, e, i, raan, argp, nu = coe.a, coe.e, coe.i, coe.raan, coe.argp, coe.nu
    p = a * (1 - e**2)
    r_pf = p / (1 + e * np.cos(nu)) * np.array([np.cos(nu), np.sin(nu), 0.0])
    v_pf = np.sqrt(mu / p) * np.array([-np.sin(nu), e + np.cos(nu), 0.0])

    def rot3(t):
        c, s = np.cos(t), np.sin(t)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

    def rot1(t):
        c, s = np.cos(t), np.sin(t)
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

    Q = rot3(raan) @ rot1(i) @ rot3(argp)
    return Q @ r_pf, Q @ v_pf


def propagate_kepler(r0, v0, dt, mu=MU_EARTH):
    """Analytic two-body propagation of an elliptic orbit by ``dt`` seconds."""
    coe = rv_to_coe(r0, v0, mu)
    e = coe.e
    n = np.sqrt(mu / coe.a**3)
    E0 = np.arctan2(np.sqrt(1 - e**2) * np.sin(coe.nu), e + np.cos(coe.nu))
    M = E0 - e * np.sin(E0) + n * dt
    E = solve_kepler(M, e)
    nu = np.arctan2(np.sqrt(1 - e**2) * np.sin(E), np.cos(E) - e)
    return coe_to_rv(OrbitalElements(coe.a, e, coe.i, coe.raan, coe.argp, nu), mu)
