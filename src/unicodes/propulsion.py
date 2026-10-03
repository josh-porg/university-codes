"""Rocket-assisted projectile description (AE 721 Missile Design)."""

from __future__ import annotations

import math
from dataclasses import dataclass, field


@dataclass
class RocketBullet:
    """Design parameters of a rocket-assisted bullet. SI units throughout."""

    name: str
    hydraulic_diameter: float  # body diameter (m)
    nozzle_exit_diameter: float  # m
    projectile_length: float  # m
    nose_length: float  # m
    nose_tip_diameter: float  # m
    boattail_diameter: float  # m (= hydraulic_diameter if no boat-tail)
    projectile_mass: float  # wet mass (kg)
    specific_impulse: float  # s
    propellant_mass: float  # kg
    motor_on_time: float  # s
    motor_off_time: float  # s
    nozzle_exit_area: float = field(init=False)  # m^2

    def __post_init__(self):
        self.nozzle_exit_area = math.pi * (self.nozzle_exit_diameter / 2) ** 2

    @property
    def burn_time(self) -> float:
        return self.motor_off_time - self.motor_on_time

    @property
    def dry_mass(self) -> float:
        return self.projectile_mass - self.propellant_mass
