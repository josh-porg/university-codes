"""AE 546 Homework 4, problems 4.6 and 4.7: thin-airfoil theory for a two-piece camber line."""

import numpy as np

from unicodes.aero.thin_airfoil import two_piece_parabolic_camber

af = two_piece_parabolic_camber(0.4, 0.25, 0.0, 0.8, 0.111, 0.2, 0.8)
alpha = np.deg2rad(3)
print(f"alpha_L0 = {np.rad2deg(af.zero_lift_angle):.4f} deg")
print(f"c_l(3 deg) = {af.c_l(alpha):.4f}")
print(f"c_m,c/4 = {af.c_m_quarter_chord:.4f}")
print(f"x_cp/c = {af.center_of_pressure(alpha):.4f}")
