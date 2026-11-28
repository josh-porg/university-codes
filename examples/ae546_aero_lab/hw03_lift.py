"""AE 546 Homework 3, problems 4.1, 4.2, 4.11: section lift per unit span."""

from unicodes.atmosphere import isa
from unicodes.units import SLUG_PER_FT3

c, V, c_l = 2.0, 50.0, 0.65
rho_bg = isa(0).rho / SLUG_PER_FT3
print("4.1: L' =", c_l * 0.5 * rho_bg * V**2 * c, "lbf/ft")
rho = isa(0).rho
print("4.2: c_l =", 1353 / (0.5 * rho * V**2 * c))
rho = isa(3000).rho
print("4.11: L' =", c_l * 0.5 * rho * 60.0**2 * c, "N/m")
