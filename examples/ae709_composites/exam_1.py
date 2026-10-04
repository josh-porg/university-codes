"""AE 709 exam 1 problems 4, 5 and 7 (``Exam_1_Problem_4``, ``Exam_1_Problem_5``, ``Exam_1_Problem_7``).

* Problem 4: micromechanics (rule of mixtures, Halpin-Tsai) at three fibre
  volume fractions and E_x of [45/90/-45/0/90]_s. The MATLAB passed
  ``Ply(E_1, E_2, nu_12, G_12)``, swapping G_12 and nu_12; fixed here.
* Problem 5: E_x, E_y, G_xy, nu_xy of [0/90/45/-45] even-symmetric,
  odd-symmetric and unsymmetric.
* Problem 7: sensitivities of the in-plane properties to the 0/45/90 ply
  fractions (finite differences instead of the symbolic derivatives) and the
  19 laminate families around 50/40/10.
"""

import numpy as np

from unicodes.composites import Laminate, Ply, mechanical_properties

# Problem 4 (msi, in)
for V_f in (0.45, 0.5, 0.55):
    E1, E2, nu12, G12 = mechanical_properties(25, 7.5, 0.3, 0.5, 0.25, 0.3, V_f, 2, 1)
    ply = Ply(E1, E2, G12, nu12)
    lam = Laminate([ply], None, 0.01, np.deg2rad([45, 90, -45, 0, 90]), 2)
    print(f"P4 V_f {V_f}: E_1 {E1:.3f}, E_2 {E2:.3f}, nu_12 {nu12:.3f}, G_12 {G12:.3f} msi; "
          f"t {lam.t_laminate:.2f} in, E_x {lam.E_x:.3f} msi")

# Problem 5
ply = Ply(18.7, 1.9, 0.85, 0.3)
for sym, label in ((2, "_s"), (1, "_os"), (0, "")):
    lam = Laminate([ply], None, 0.01, np.deg2rad([0, 90, 45, -45]), sym)
    print(f"P5 [0/90/45/-45]{label}: E_x {lam.E_x:.3f}, E_y {lam.E_y:.3f}, G_xy {lam.G_xy:.3f}, nu_xy {lam.nu_xy:.4f}")


# Problem 7: membrane properties of a 0/+-45/90 family with thickness fractions t0, t45, t90
def props(t0, t45, t90):
    A = ply.Qbar(0) * t0 + (ply.Qbar(np.pi / 4) + ply.Qbar(-np.pi / 4)) / 2 * t45 + ply.Qbar(np.pi / 2) * t90
    a = np.linalg.inv(A)
    t = t0 + t45 + t90
    return np.array([1 / (a[0, 0] * t), 1 / (a[1, 1] * t), 1 / (a[2, 2] * t), -a[0, 1] / a[0, 0]])


base, h = np.array([0.5, 0.4, 0.1]), 1e-6
sens = np.column_stack([(props(*(base + h * e)) - props(*(base - h * e))) / (2 * h) for e in np.eye(3)])
print("P7 sensitivities d[E_x, E_y, G_xy, nu_xy]/d[t0, t45, t90] at 50/40/10:\n", np.round(sens, 3))
families = []
for i in range(1, 6):
    R0 = 48 + (2 * i - 4)
    for R45 in range(max(36, 42 - 2 * i), min(44, 50 - 2 * i) + 1, 2):
        families.append((R0, R45, 100 - R0 - R45))
res = np.array([[Laminate([ply], None, [r0, r45 / 2, r45 / 2, r90], np.deg2rad([0, 45, -45, 90]), 2).__dict__[k]
                 for k in ("E_x", "E_y", "G_xy", "nu_xy")] for r0, r45, r90 in families])
nominal = res[families.index((50, 40, 10))]
print("P7 family (0/45/90 %)   E_x     E_y     G_xy    nu_xy   max |change|")
for fam, r in zip(families, res):
    print(f"   {fam!s:16s} {r[0]:7.3f} {r[1]:7.3f} {r[2]:7.3f} {r[3]:7.4f}  {np.max(np.abs((r - nominal) / nominal)):.3f}")
