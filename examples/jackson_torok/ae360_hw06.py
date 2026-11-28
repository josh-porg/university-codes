"""Jackson Torok, AE 360 homework 6 (``AE360_Torok_HW6``): Kepler's equation by Newton iteration.

Eccentric anomaly for mean anomalies 0-350 deg and e = 0.01, 0.1, 0.5,
0.9, starting from ``E0 = M + e`` (``M - e`` for M in (-pi, 0) or M > pi).
The MATLAB stopped when the change of the step fell below the tolerance;
here it stops when the step itself does, and the result is checked
against :func:`unicodes.orbital.solve_kepler`.
"""

import numpy as np

from unicodes.orbital import solve_kepler


def kepler_newton(M, e, tol=1e-6):
    E = M - e if (-np.pi < M < 0) or M > np.pi else M + e
    while True:
        dE = (M - E + e * np.sin(E)) / (1 - e * np.cos(E))
        E += dE
        if abs(dE) < tol:
            return E


print(" M (deg)     e     E (deg)")
for e in (0.01, 0.1, 0.5, 0.9):
    for M_deg in (0, 30, 90, 135, 180, 200, 350):
        M = np.radians(M_deg)
        E = kepler_newton(M, e)
        assert abs(np.mod(E - solve_kepler(M, e) + np.pi, 2 * np.pi) - np.pi) < 1e-6
        print(f"{M_deg:7d} {e:6.2f} {np.degrees(E):10.4f}")
