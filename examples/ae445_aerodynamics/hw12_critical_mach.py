"""AE 445 Homework 12, problems 3.3 and 3.7: critical Mach number of an airfoil.

Reads the upper-surface C_p data the course supplied
(``Prob_3_1_cp_upper.xlsx``; pass its path) and finds where the
Prandtl-Glauert and Karman-Tsien corrected minimum C_p meets the sonic line.
"""

import sys

import numpy as np

from unicodes.aero import wing
from unicodes.atmosphere import speed_of_sound
from unicodes.gasdynamics import critical_pressure_coefficient
from unicodes.units import MPH


def main(cp_upper):
    cp0 = float(np.min(cp_upper))
    for method in ("prandtl_glauert", "karman_tsien"):
        print(f"M_crit ({method}) = {wing.critical_mach_from_cp(cp0, method):.4f}")

    M = np.linspace(0.01, 0.99, 500)
    try:
        import matplotlib.pyplot as plt

        plt.plot(M, critical_pressure_coefficient(M), label="Sonic line")
        plt.plot(M, wing.prandtl_glauert(cp0, M), label="Prandtl-Glauert")
        plt.plot(M, wing.karman_tsien(cp0, M), label="Karman-Tsien")
        plt.ylim(0, -1.5)
        plt.xlabel(r"$M_\infty$")
        plt.ylabel("$C_p$")
        plt.legend()
        plt.show()
    except ImportError:
        pass

    # 3.7: 450 mph at sea level
    M = 450 * MPH / float(speed_of_sound(288.2))
    print(f"3.7: M = {M:.3f}, Cp_PG = {wing.prandtl_glauert(cp0, M):.4f}, Cp_KT = {wing.karman_tsien(cp0, M):.4f}")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        import pandas as pd  # reading .xlsx needs pandas + openpyxl

        main(pd.read_excel(sys.argv[1], header=None).to_numpy().ravel())
    else:
        print("usage: python hw12_critical_mach.py Prob_3_1_cp_upper.xlsx")
