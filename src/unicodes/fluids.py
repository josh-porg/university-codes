"""Incompressible pipe-flow helpers (AE 345)."""

from __future__ import annotations

from dataclasses import dataclass

G0 = 9.80665  # m/s^2


@dataclass
class PipeSection:
    """A pipe fitting with a minor-loss coefficient between two areas.

    Give either ``V_initial`` or ``V_final``; the other follows from
    incompressible continuity ``A1 V1 = A2 V2``.
    """

    area_initial: float  # m^2
    area_final: float  # m^2
    K_L: float  # minor-loss coefficient (based on the inlet velocity)
    rho: float  # kg/m^3
    V_initial: float | None = None  # m/s
    V_final: float | None = None  # m/s

    def __post_init__(self):
        if self.V_initial is None and self.V_final is None:
            raise ValueError("Either V_initial or V_final must be provided")
        if self.V_final is None:
            self.V_final = self.area_initial * self.V_initial / self.area_final
        elif self.V_initial is None:
            self.V_initial = self.area_final * self.V_final / self.area_initial

    @property
    def pressure_loss(self) -> float:
        """Minor pressure loss K_L * rho * V1^2 / 2 (Pa)."""
        return self.K_L * self.rho * self.V_initial**2 / 2

    def head_loss(self, g: float = G0) -> float:
        """Minor head loss K_L * V1^2 / (2 g) (m)."""
        return self.K_L * self.V_initial**2 / (2 * g)

    @property
    def delta_p(self) -> float:
        """Static pressure change p2 - p1 (Pa) from Bernoulli plus the minor loss."""
        return self.rho * (self.V_initial**2 - self.V_final**2) / 2 - self.pressure_loss


def reynolds_number(rho, V, L, mu):
    return rho * V * L / mu


def friction_factor_laminar(Re):
    """Darcy friction factor of fully developed laminar pipe flow, ``64 / Re``."""
    return 64 / Re


def friction_factor_colebrook(Re, rel_roughness):
    """Darcy friction factor from the Colebrook equation (``getfrictionfactorCoelbrook``)."""
    from scipy.optimize import brentq

    def g(inv_sqrt_f):
        return inv_sqrt_f + 2 * __import__("math").log10(rel_roughness / 3.7 + 2.51 * inv_sqrt_f / Re)

    return 1 / brentq(g, 0.5, 100) ** 2


def friction_factor_haaland(Re, rel_roughness):
    """Haaland's explicit approximation (``getfrictionfactorHaaland``; ``CoelbrookAlt`` is the same formula)."""
    import math

    return (1 / (-1.8 * math.log10((rel_roughness / 3.7) ** 1.11 + 6.9 / Re))) ** 2


def darcy_pressure_drop(f, length, diameter, rho, V):
    """``f (L / D) rho V^2 / 2``."""
    return f * length / diameter * rho * V**2 / 2
