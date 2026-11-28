"""Blade-element momentum (BEMT) analysis of a propeller (AE 722).

The blade has linear taper and the twist distribution
``theta(x) = theta_0 + (theta_1 + theta_tw * x) * x**twist_exponent``
with ``x = r / R``. Inflow uses the Prandtl tip-loss factor. Lift is
linear up to stall; drag is ``C_D0 + d1 alpha + d2 alpha^2``.

    blade = PropellerBlade(np.deg2rad(45), np.deg2rad(15), np.deg2rad(-42), -0.7, 0.1, radius=1.1)
    perf = blade.analyze(omega=134.5, V=75 * 0.514444, rho=1.099)
    perf = blade.match_power(19.7e3, V=..., rho=...)   # finds the RPM
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq


@dataclass(frozen=True)
class PropellerPerformance:
    omega: float  # rad/s
    thrust: float  # N
    shaft_power: float  # W
    aero_power: float  # W, thrust x flight speed
    efficiency: float
    C_T: float
    C_P: float
    alpha: np.ndarray  # section angle of attack along the span (rad), stall-limited


@dataclass(frozen=True)
class PropellerBlade:
    theta_0: float  # rad
    theta_1: float  # rad
    theta_tw: float  # rad
    twist_exponent: float
    root_chord: float  # m
    radius: float  # m
    n_blades: int = 3
    taper: float = 0.8  # tip chord / root chord
    root_cutout: float = 0.15  # fraction of radius
    c_l_alpha: float = 0.108 * 180 / np.pi  # NACA 0012 (1/rad)
    C_D0: float = 0.01
    d1: float = 0.025
    d2: float = 0.65
    alpha_stall: float = np.deg2rad(17)  # stall angle above zero lift (rad)
    n_segments: int = 100
    tip_loss_iterations: int = 4

    def twist(self, x):
        return self.theta_0 + (self.theta_1 + self.theta_tw * x) * x**self.twist_exponent

    def _sections(self, V, V_tip):
        dx = 1 / self.n_segments
        x_out = np.arange(1, self.n_segments + 1) * dx
        x = x_out - dx / 2  # mid-span of each segment
        chord = self.root_chord * (1 - x_out * (1 - self.taper))
        sigma = np.where(
            x > self.root_cutout, self.n_blades * chord / (2 * np.pi * x_out * self.radius), 1e-4
        )
        lam_c = V / V_tip
        theta = self.twist(x)
        a = self.c_l_alpha

        def inflow(F):
            b = sigma * a / (16 * F) - lam_c / 2
            return np.sqrt(b**2 + sigma * a * theta * x / (8 * F)) - b

        lam = inflow(1.0)
        for _ in range(self.tip_loss_iterations):
            F = 2 / np.pi * np.arccos(np.exp(-self.n_blades / 2 * (1 - x) / lam))
            lam = inflow(F)

        alpha = np.minimum(theta - np.arctan2(lam, x), self.alpha_stall)
        C_d = self.C_D0 + self.d1 * alpha + self.d2 * alpha**2
        dC_T = sigma * a / 2 * (theta * x**2 - lam * x) * dx
        if np.degrees(alpha[0]) >= np.degrees(self.alpha_stall):
            dC_T[0] = sigma[0] * a / 2 * self.alpha_stall * dx
        dC_P = lam * dC_T + sigma / 2 * C_d * x**3 * dx
        return dC_T.sum(), dC_P.sum(), alpha

    def analyze(self, omega: float, V: float, rho: float) -> PropellerPerformance:
        """Performance at rotation rate ``omega`` (rad/s), flight speed ``V`` (m/s), density ``rho``."""
        V_tip = omega * self.radius
        C_T, C_P, alpha = self._sections(V, V_tip)
        area = np.pi * self.radius**2
        T = C_T * rho * area * V_tip**2
        P_shaft = C_P * rho * area * V_tip**3
        P_aero = T * V
        return PropellerPerformance(omega, T, P_shaft, P_aero, P_aero / P_shaft, C_T, C_P, alpha)

    def match_power(
        self, aero_power: float, V: float, rho: float, omega_bounds=(20.0, 1000.0)
    ) -> PropellerPerformance:
        """Find the rotation rate that delivers ``aero_power`` (thrust power, W)."""
        omega = brentq(lambda w: self.analyze(w, V, rho).aero_power - aero_power, *omega_bounds, xtol=1e-6)
        return self.analyze(omega, V, rho)
