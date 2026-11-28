"""AE 550 homework 2-4: rotations, flight-path and air-flow angles with wind (``HW_2``, ``HW_3``).

HW 4 (symbolic linearisation of the rigid-body equations) is a derivation;
its result is the small-perturbation model used by
:func:`unicodes.flight_dynamics.sixdof.linearize`.
"""

import numpy as np

from unicodes.flight_dynamics import kinematics as k

d = np.rad2deg
# HW 2 problem 4
att = np.deg2rad([10, 5, -45])
V_B = np.array([850.0, -30, 20])
V_I = k.dcm_inertial_to_body(att).T @ V_B
g, chi = k.flight_path_angles(V_I)
print(f"HW2 P4: V_I = {V_I}, gamma {d(g):.3f}, heading {d(chi):.3f} deg")
# problems 5 and 6
V_I = k.velocity_from_path_angles(450, np.deg2rad(50), np.deg2rad(25))
vB1 = k.dcm_inertial_to_body(np.deg2rad([2, 60, 28])) @ V_I
vB2 = k.dcm_inertial_to_body(np.deg2rad([4, 55, 33])) @ V_I
print(f"HW2 P5/6: V_I = {V_I}, V_B = {vB1}, guessed attitude error = {vB1 - vB2}")
# HW 3
V_I = np.array([850.0, -30, 20])
H = k.dcm_inertial_to_body(np.deg2rad([20, 5, -5]))
for wind in ([0, 0, 0], [-80 * np.sin(np.pi / 4), 80 * np.sin(np.pi / 4), 0], [80, 0, 0]):
    a, b = k.air_flow_angles(H @ (V_I - np.array(wind)))
    print(f"HW3 wind {np.round(wind, 1)}: alpha {d(a):.3f}, beta {d(b):.3f} deg")
