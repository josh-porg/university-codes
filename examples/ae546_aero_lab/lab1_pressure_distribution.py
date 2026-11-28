"""AE 546 Lab 1: lift, drag and moment of a NACA 23012 from pressure taps.

Expects the lab spreadsheet ``Lab1_Data.xlsx`` (A1:G51): row 1 angle of
attack (deg), row 2 static head, row 3 total head, then 48 tap heads, upper
surface from the trailing edge first. Replaces ``Lab_1``,
``Lab_1_preliminary``, ``Lab_1_version_two...`` and ``lab1_ploter_maybe``.
"""

import sys

import numpy as np

from unicodes.aero import naca, pressure

MU, RHO_AIR, RHO_LIQUID, CHORD = 1.78e-5, 1.225, 836.0, 0.195


def analyse(data):
    alpha = np.deg2rad(data[0])
    h_s, h_t, h = data[1], data[2], data[3:]
    V = pressure.manometer_velocity(h_t, h_s, RHO_LIQUID, RHO_AIR)
    Re = RHO_AIR * V * CHORD / MU
    cp = pressure.cp_from_heads(h, h_s, h_t)
    af = naca.naca5("23012", n=24)
    x_u, z_u = af.upper
    upper = (x_u[::-1], z_u[::-1])  # leading edge -> trailing edge
    cp_upper = cp[:24][::-1]  # taps were numbered from the trailing edge
    cp_lower = cp[24:]
    return Re, cp, pressure.section_coefficients(alpha, cp_upper, cp_lower, upper, af.lower)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: python lab1_pressure_distribution.py Lab1_Data.xlsx")
    import pandas as pd

    data = pd.read_excel(sys.argv[1], header=None).to_numpy(dtype=float)[:51, :7]
    Re, cp, res = analyse(data)
    for a, r, cl, cd, cm in zip(data[0], Re, res.c_l, res.c_d, res.c_m_qc):
        print(f"alpha {a:5.1f}  Re {r:9.0f}  c_l {cl:7.3f}  c_d {cd:7.4f}  c_m,c/4 {cm:7.3f}")
