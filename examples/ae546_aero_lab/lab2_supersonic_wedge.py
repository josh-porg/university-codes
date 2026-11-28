"""AE 546 Lab 2: oblique shock then expansion over a 5 deg wedge (``Lab_2_driver``, ``Lab_2_function``).

Compares the free-stream Mach number inferred from the measured shock angle
with the one from the nozzle area ratio, as in the lab report tables.
"""

import numpy as np

from unicodes import gasdynamics as gd
from unicodes.io import to_latex_table

R, P1, T1, THETA = 287.052874, 101325.0, 288.0, np.deg2rad(5)


def regions(M1, beta):
    """p, p0, T, T0, Mach angle and Mach in the free stream, behind the shock and after the expansion."""
    iso1 = gd.isentropic(M1)
    p01, T01 = P1 * iso1.p0_p, T1 * iso1.T0_T
    s = gd.oblique_shock(M1, beta, THETA)
    M2, p2, T2 = float(s.M2), P1 * s.p2_p1, T1 * s.T2_T1
    p02 = p01 * s.p02_p01
    M3 = gd.expansion_fan(M2, THETA)
    iso3 = gd.isentropic(M3)
    p3, T3 = p02 / iso3.p0_p, T01 / iso3.T0_T
    M = np.array([M1, M2, M3])
    return np.column_stack([[P1, p2, p3], [p01, p02, p02], [T1, T2, T3], [T01, T01, T01], gd.mach_angle(M), M])


M_from_beta = 2.0  # from the measured 30 deg shock angle
M_from_area = gd.mach_from_area_ratio(3.4 * 2 / (0.875 * 2), supersonic=True)
headers = ["p", "p_0", "T", "T_0", r"\mu", "M"]
A = regions(M_from_beta, np.deg2rad(30))
B = regions(M_from_area, gd.wave_angle(M_from_area, THETA))
names = ["1", "2", "3"]
print(f"M from shock angle {M_from_beta}, from area ratio {M_from_area:.4f}")
print(to_latex_table(A, headers, row_names=names))
print(to_latex_table(B, headers, row_names=names))
print(to_latex_table(np.abs((A - B) / A) * 100, headers, row_names=names))
