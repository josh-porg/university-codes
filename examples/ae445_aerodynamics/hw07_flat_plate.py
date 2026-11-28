"""AE 445 Homework 7, problem 2.8: drag of a 3 ft x 6 ft plate at 100 mph, sea level.

Laminar up to Re = 750,000, turbulent after; both sides wetted.
"""

from unicodes.aero import boundary_layer as bl
from unicodes.atmosphere import isa
from unicodes.units import FT, LBF, MPH

chord, span, V = 3 * FT, 6 * FT, 100 * MPH
air = isa(0)
D_laminar_all = bl.flat_plate_drag(V, chord, span, air.rho, air.nu, Re_crit=1e12)
D_mixed = bl.flat_plate_drag(V, chord, span, air.rho, air.nu, Re_crit=750_000)
D_turbulent = bl.flat_plate_drag(V, chord, span, air.rho, air.nu)
x_t = bl.transition_distance(750_000, V, air.nu)
print(f"transition at x = {x_t / FT:.3f} ft")
print(f"all laminar   D = {D_laminar_all / LBF:.4f} lbf")
print(f"mixed         D = {D_mixed / LBF:.4f} lbf")
print(f"all turbulent D = {D_turbulent / LBF:.4f} lbf")
