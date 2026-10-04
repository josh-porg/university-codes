"""Kepler's equation by Newton's method for a table of M and e (``Joshua_Poznanski_HW6``, ``newtons_method_experimentation``)."""

import numpy as np

from unicodes.orbital import solve_kepler

M = np.deg2rad([0, 30, 90, 135, 180, 200, 350])
for e in (0.01, 0.1, 0.5, 0.9):
    E = [solve_kepler(m, e) for m in M]
    print(f"e = {e}: E (deg) =", np.round(np.rad2deg(np.mod(E, 2 * np.pi)), 4))
