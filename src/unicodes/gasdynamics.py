"""Compressible-flow relations for a calorically perfect gas (AE 546, AE 573).

Isentropic ratios, normal and oblique shocks, Prandtl-Meyer expansions,
pitot-tube Mach number and inlet figures of merit. All functions take
``gamma`` (default 1.4) and work elementwise on arrays. Angles in radians.

The MATLAB versions solved the Prandtl-Meyer inverse and the supersonic
pitot formula symbolically (``vpasolve``/``solve``); here they use a
bracketed root finder.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

GAMMA = 1.4


@dataclass
class IsentropicRatios:
    """Total-to-static ratios and the area ratio A/A* at a Mach number."""

    p0_p: np.ndarray
    rho0_rho: np.ndarray
    T0_T: np.ndarray
    A_Astar: np.ndarray


def isentropic(M, gamma=GAMMA) -> IsentropicRatios:
    """Isentropic total/static ratios and A/A* at Mach ``M`` (``isentropicFlowProperties``)."""
    M = np.asarray(M, dtype=float)
    t = 1 + (gamma - 1) / 2 * M**2
    with np.errstate(divide="ignore"):
        A = np.sqrt(1 / M**2 * (2 / (gamma + 1) * t) ** ((gamma + 1) / (gamma - 1)))
    return IsentropicRatios(t ** (gamma / (gamma - 1)), t ** (1 / (gamma - 1)), t, A)


def total_pressure(p, M, gamma=GAMMA):
    """Total pressure from static pressure and Mach (``totalPressureFromStaticAndMach``)."""
    return p * isentropic(M, gamma).p0_p


def mach_from_area_ratio(A_Astar, supersonic: bool, gamma=GAMMA):
    """Invert A/A* for Mach on the subsonic or supersonic branch."""
    if A_Astar < 1:
        raise ValueError("A/A* must be >= 1")
    if np.isclose(A_Astar, 1):
        return 1.0
    f = lambda M: isentropic(M, gamma).A_Astar - A_Astar  # noqa: E731
    return brentq(f, 1 + 1e-12, 100) if supersonic else brentq(f, 1e-9, 1 - 1e-12)


@dataclass
class ShockResult:
    """Downstream state across a shock, as ratios downstream/upstream."""

    M2: np.ndarray
    rho2_rho1: np.ndarray
    p2_p1: np.ndarray
    T2_T1: np.ndarray
    p02_p01: np.ndarray  # total pressure recovery


def normal_shock(M1, gamma=GAMMA) -> ShockResult:
    """Normal-shock jump relations for upstream Mach ``M1`` > 1 (``normalShockProperties``)."""
    M1 = np.asarray(M1, dtype=float)
    M2 = np.sqrt(((gamma - 1) * M1**2 + 2) / (2 * gamma * M1**2 - (gamma - 1)))
    rho = (gamma + 1) * M1**2 / (2 + (gamma - 1) * M1**2)
    p = 1 + 2 * gamma / (gamma + 1) * (M1**2 - 1)
    p0 = rho ** (gamma / (gamma - 1)) * ((gamma + 1) / (2 * gamma * M1**2 - (gamma - 1))) ** (
        1 / (gamma - 1)
    )
    return ShockResult(M2, rho, p, p / rho, p0)


def pitot_rayleigh(M1, gamma=GAMMA):
    """Pitot pressure over upstream static, p02/p1, behind a normal shock."""
    M1 = np.asarray(M1, dtype=float)
    return normal_shock(M1, gamma).p02_p01 * isentropic(M1, gamma).p0_p


def oblique_shock(M1, beta, theta=None, gamma=GAMMA) -> ShockResult:
    """Oblique-shock relations for wave angle ``beta`` (``obliqueShockProperties``).

    ``theta`` (deflection) is computed from theta-beta-M if not given.
    """
    M1 = np.asarray(M1, dtype=float)
    if theta is None:
        theta = deflection_angle(M1, beta, gamma)
    n = normal_shock(M1 * np.sin(beta), gamma)
    return ShockResult(n.M2 / np.sin(beta - theta), n.rho2_rho1, n.p2_p1, n.T2_T1, n.p02_p01)


def deflection_angle(M1, beta, gamma=GAMMA):
    """Flow deflection theta for wave angle ``beta`` (theta-beta-M relation)."""
    M1 = np.asarray(M1, dtype=float)
    return np.arctan(
        2 / np.tan(beta) * (M1**2 * np.sin(beta) ** 2 - 1) / (M1**2 * (gamma + np.cos(2 * beta)) + 2)
    )


def wave_angle(M1, theta, gamma=GAMMA, strong: bool = False):
    """Shock wave angle beta for deflection ``theta``; NaN if the shock detaches."""
    mu = mach_angle(M1)
    betas = np.linspace(mu, np.pi / 2, 2001)
    th = deflection_angle(M1, betas, gamma)
    k = int(np.argmax(th))
    if theta > th[k]:
        return np.nan
    if theta == 0:
        return np.pi / 2 if strong else mu
    lo, hi = (betas[k], np.pi / 2) if strong else (mu, betas[k])
    return brentq(lambda b: deflection_angle(M1, b, gamma) - theta, lo, hi)


def mach_angle(M):
    """Mach wave angle asin(1/M) (``machWaveAngle``)."""
    return np.arcsin(1 / np.asarray(M, dtype=float))


def prandtl_meyer(M, gamma=GAMMA):
    """Prandtl-Meyer angle nu(M) for M >= 1 (``prandtlMeyerFunction``)."""
    M = np.asarray(M, dtype=float)
    g = np.sqrt((gamma + 1) / (gamma - 1))
    return g * np.arctan(np.sqrt((M**2 - 1) / g**2)) - np.arctan(np.sqrt(M**2 - 1))


def max_turning_angle(gamma=GAMMA):
    """Maximum Prandtl-Meyer angle (expansion to vacuum)."""
    return np.pi / 2 * (np.sqrt((gamma + 1) / (gamma - 1)) - 1)


def inverse_prandtl_meyer(nu, gamma=GAMMA):
    """Mach number whose Prandtl-Meyer angle is ``nu`` (rad)."""
    if not 0 <= nu < max_turning_angle(gamma):
        raise ValueError("nu outside the Prandtl-Meyer range")
    if nu == 0:
        return 1.0
    return brentq(lambda M: prandtl_meyer(M, gamma) - nu, 1.0, 1e4, xtol=1e-14)


def expansion_fan(M1, theta, gamma=GAMMA):
    """Mach number after turning supersonic flow through ``theta`` (``expansionFanProperties``)."""
    return inverse_prandtl_meyer(prandtl_meyer(M1, gamma) + theta, gamma)


def mach_from_pitot_subsonic(p0, p, gamma=GAMMA):
    """Mach from pitot total and static pressure in subsonic flow (``MachFromPitoCompresible``)."""
    return np.sqrt(2 / (gamma - 1) * ((np.asarray(p0) / p) ** ((gamma - 1) / gamma) - 1))


def mach_from_pitot_supersonic(p02, p1, gamma=GAMMA):
    """Upstream Mach from pitot pressure behind a normal shock (Rayleigh pitot formula).

    ``p02`` is the pitot (total) pressure, ``p1`` the free-stream static pressure
    (``MachFromPitoSupersonic``).
    """
    ratio = p02 / p1
    if ratio < pitot_rayleigh(1.0, gamma):
        raise ValueError("pitot ratio too low for supersonic flow")
    return brentq(lambda M: pitot_rayleigh(M, gamma) - ratio, 1.0, 100.0)


def pressure_coefficient(p, p_inf, p0):
    """C_p = (p - p_inf)/(p0 - p_inf) (``CpFromPressure``; same form as ``CpFromHead``)."""
    return (np.asarray(p) - p_inf) / (p0 - p_inf)


def critical_pressure_coefficient(M, gamma=GAMMA):
    """C_p at which the local flow first reaches Mach 1, for free-stream ``M`` (``getCpCrit``).

    The AE 445 HW 12 script passed the specific *weight* of air as ``gamma``;
    it is the ratio of specific heats.
    """
    M = np.asarray(M, dtype=float)
    return 2 / (gamma * M**2) * ((((gamma - 1) * M**2 + 2) / (gamma + 1)) ** (gamma / (gamma - 1)) - 1)


# Inlets (AE 573)


def inlet_total_pressure_recovery(p0_freestream, p0_engine_face):
    """pi_d = pt2/pt0."""
    return p0_engine_face / p0_freestream


def inlet_adiabatic_efficiency(M0, p0, pt2, gamma=GAMMA):
    """Inlet adiabatic (kinetic energy) efficiency eta_d."""
    return ((pt2 / p0) ** ((gamma - 1) / gamma) - 1) / (M0**2 * (gamma - 1) / 2)


def inlet_entropy_rise(pi_d):
    """Non-dimensional entropy rise Delta s / R = -ln(pi_d)."""
    return -np.log(pi_d)
