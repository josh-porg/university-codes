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
