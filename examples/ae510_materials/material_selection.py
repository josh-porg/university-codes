"""AE 510 materials and processes: HW 4 member weights and the supplemental-project beam material
screen (``HW_4``, ``save_materials``, ``supplimental_project_one*``).

The project screened every material of ``MaterialDataSP1`` (an Excel export
holding E [MPa], F_TU [ksi], F_TY [ksi], rho [lbf/in^3] and cost [$/lbf];
not in the repository) for a 6 in x 1 in x 0.25 in beam carrying p = 500 lbf::

    python material_selection.py MaterialDataSP1.csv

Without a file a small built-in table is screened. Failure checks follow the
MATLAB: bending stress 120 p against F_TU and 1.5 F_TY, deflection 2000 p / E
against 0.08 in, and cost against $1. (``supplimental_project_one-ENGR-*``
are unfinished drafts of the same script.)
"""

import argparse
import csv

import numpy as np

from unicodes.materials import Material

# HW 4: mass and weight of a 5000 N tension member in six materials
sigma = np.array([20, 51, 41, 40, 41, 77]) * 1e6
rho = np.array([0.45, 1.7, 1.4, 1.1, 0.98, 1.2]) * 1e3
A = 5000 / sigma
print("HW 4: area (m^2)", np.round(A, 8), "\n      weight per metre (N/m)", np.round(rho * A * 9.81, 4))

p = argparse.ArgumentParser()
p.add_argument("csv", nargs="?", help="columns Material, EMPa, F_TUksi, F_TYksi, rholbfin3, lbf")
a = p.parse_args()

# Imperial units as in the project (psi, lbf/in^3, $/lbf)
if a.csv:
    with open(a.csv, newline="") as f:
        materials = [Material(r["Material"], float(r["EMPa"]) * 1e6 / 6894.76, float(r["F_TUksi"]) * 1e3,
                              float(r["F_TYksi"]) * 1e3, float(r["rholbfin3"]), float(r["lbf"])) for r in csv.DictReader(f)]
else:
    materials = [Material("Al 6061-T6", 10e6, 45e3, 40e3, 0.098, 2.0), Material("Steel 4130", 29.7e6, 97e3, 63e3, 0.284, 1.5),
                 Material("Douglas fir", 1.9e6, 12e3, 7e3, 0.019, 1.0), Material("Ti-6Al-4V", 16.5e6, 138e3, 128e3, 0.16, 25.0)]

volume = 0.25 * 1 * 6  # in^3 (the MATLAB multiplied by (1/12)^3, converting to ft^3 while rho is per in^3)
load = 500.0
sigma_bend = 120 * load
for m in materials:
    weight = volume * m.rho
    cost = m.cost_density * weight
    if sigma_bend > m.F_tu:
        verdict = "fails: ultimate bending stress"
    elif sigma_bend > 1.5 * m.F_ty:
        verdict = "fails: yield bending stress"
    elif 2000 * load / m.E > 0.08:
        verdict = "fails: deflection"
    elif cost > 1:
        verdict = "fails: cost"
    else:
        verdict = "passes"
    print(f"{m.name:14s} weight {weight:6.3f} lbf, cost ${cost:6.2f}: {verdict}")
