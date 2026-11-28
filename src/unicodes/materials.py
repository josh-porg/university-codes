"""Basic structural material record (AE 510 Materials and Processes)."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Material:
    """Isotropic material properties for quick sizing and trade studies."""

    name: str
    E: float  # Young's modulus (Pa)
    F_tu: float  # ultimate tensile strength (Pa)
    F_ty: float  # tensile yield strength (Pa)
    rho: float  # density (kg/m^3)
    cost_density: float  # cost per unit mass (e.g. $/kg)

    @property
    def specific_stiffness(self) -> float:
        return self.E / self.rho

    @property
    def specific_strength(self) -> float:
        return self.F_tu / self.rho
