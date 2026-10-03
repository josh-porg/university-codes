"""AE 445 Homework 4, problem 6.46: free-vortex (hurricane) velocities."""

import numpy as np

from unicodes.aero.wing import free_vortex_velocity

v_max, r_wall = 160.0, 10.0  # mph, mi
circulation = -2 * np.pi * v_max * r_wall  # mi^2/hr
r = np.array([10, 20, 30, 40])  # mi
print("Gamma =", circulation, "mi^2/hr")
print("V(r) =", free_vortex_velocity(r, circulation), "mph")
