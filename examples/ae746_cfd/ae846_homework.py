"""AE 846 advanced CFD: Euler flux-Jacobian eigensystem, kappa-MUSCL cell averages and FR coefficients.

* ``HW_1_question_1``: Jacobian of the 1-D Euler flux and its eigenvalues (sympy).
* ``HW_1_question_4``: cell averages of the kappa-family quadratic reconstruction.
* ``FR_coefs``: flux-reconstruction differentiation, interpolation and correction coefficients.

``Reacting_flow_1D_equations_derivation`` is an unfinished symbolic
derivation (it ends in an undefined substitution) and is not ported.
"""

import numpy as np
import sympy as sp

from unicodes.cfd.flux_reconstruction import fr_operators

g, Q1, Q2, Q3 = sp.symbols("gamma Q_1 Q_2 Q_3")
p = (g - 1) * (Q3 - Q2**2 / (2 * Q1))
F = sp.Matrix([Q2, Q2**2 / Q1 + p, Q2 / Q1 * (Q3 + p)])
A = F.jacobian([Q1, Q2, Q3])
rho, u, pr = sp.symbols("rho u p", positive=True)
prim = {Q1: rho, Q2: rho * u, Q3: pr / (g - 1) + rho * u**2 / 2}
print("Flux Jacobian:", sp.simplify(A.subs(prim)))
print("Eigenvalues:", [sp.simplify(ev) for ev in A.subs(prim).eigenvals()])

x, xi, dx, kappa, um, ui, up = sp.symbols("x x_i Delta_x kappa u_im1 u_i u_ip1")
recon = ui + (x - xi) * (up - um) / (2 * dx) + sp.Rational(3, 2) * kappa * ((x - xi) ** 2 - dx**2 / 12) * (up - 2 * ui + um) / dx**2
for name, a, b in (("cell i", xi - dx / 2, xi + dx / 2), ("cell i-1", xi - 3 * dx / 2, xi - dx / 2)):
    print(f"Average over {name}:", sp.simplify(sp.integrate(recon, (x, a, b)) / dx))

for k in range(5):
    op = fr_operators(k)
    print(f"\nk = {k}: nodes {np.round(op.nodes, 6)}")
    print("  D =", np.round(op.D, 6).tolist())
    print("  interpolation to -1, 1 =", np.round(op.interp, 6).tolist())
    print("  dg_L at nodes =", np.round(op.dg_left, 6), " dg_R =", np.round(op.dg_right, 6))
