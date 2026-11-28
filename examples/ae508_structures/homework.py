"""AE 508 homework 6, 7 and 9: stiffness-method truss and beam problems, shear and moment diagrams.

Needs sympy (``pip install unicodes[symbolic]``) for the symbolic stiffness
solutions, as in the MATLAB Symbolic Toolbox originals.
"""

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

# HW 6: three-node truss with node 3 displaced 0.001 L along 30 deg
A, E, L = sp.symbols("A E L", positive=True)
r2 = sp.sqrt(2)
k = A * E / (2 * r2 * L) * sp.Matrix([
    [1 + r2, 1, -1, -1, -r2, 0],
    [1, 1, -1, -1, 0, 0],
    [-1, -1, 2, 0, -1, 1],
    [-1, -1, 0, 2, 1, -1],
    [-r2, 0, -1, 1, 1 + r2, -1],
    [0, 0, 1, -1, -1, 1],
])  # fmt: skip
u2 = sp.Symbol("u2")
q = sp.Matrix([0, 0, u2, 0, sp.Rational(1, 1000) * L * sp.sqrt(3) / 2, sp.Rational(1, 1000) * L / 2])
Q = k * q
u2_sol = sp.solve(Q[2], u2)[0]  # X2 = 0 (node 2 free in x)
q = q.subs(u2, u2_sol)
Q = sp.simplify(k * q)
print("HW 6: u2 =", sp.simplify(u2_sol))
print("      nodal forces Q =", list(Q))
print("      strain energy q^T k q / 2 =", sp.simplify((q.T * k * q)[0] / 2))

# HW 7 problem 1: two-element beam, v2 = 0.001 l imposed, M_z2 = 0
l = sp.Symbol("l", positive=True)
k7 = sp.Matrix([
    [3, 3 * l, -3, 3 * l, 0, 0],
    [3 * l, 4 * l**2, -3 * l, 2 * l**2, 0, 0],
    [-3, -3 * l, 9, 0, -6, 3 * l],
    # The MATLAB row 4 had 2*l and 4*l + 2*l^2 (not symmetric with column 4); corrected to the l^2 terms.
    [3 * l, 2 * l**2, 0, 6 * l**2, -3 * l, l**2],
    [0, 0, -6, -3 * l, 6, -3 * l],
    [0, 0, 3 * l, l**2, -3 * l, 2 * l**2],
])  # fmt: skip
th2 = sp.Symbol("theta_z2")
q7 = sp.Matrix([0, 0, l / 1000, th2, 0, 0])
th2_sol = sp.solve((k7 * q7)[3], th2)[0]
print("HW 7: theta_z2 =", th2_sol, " forces =", list(sp.simplify(k7 * q7.subs(th2, th2_sol))))

# HW 7 part 2: propped beam with equivalent nodal loads
p, EI = sp.symbols("p EI", positive=True)
v1, t1, Yr2 = sp.symbols("v1 theta_z1 Y_r2")
kff = EI * 32 / L**3 * sp.Matrix([[6, sp.Rational(3, 4) * L, -6], [sp.Rational(3, 4) * L, L**2 / 8, -sp.Rational(3, 4) * L],
                                  [-6, sp.Rational(3, 4) * L, sp.Rational(56, 9)]])
# The MATLAB third equivalent load lacked the factor p.
Qf = sp.Matrix([-sp.Rational(5, 32) * p, sp.Rational(9, 64) * p * L, Yr2 - sp.Rational(22, 32) * p])
sol = sp.solve(list(Qf - kff * sp.Matrix([v1, t1, 0])), [v1, t1, Yr2], dict=True)[0]
print("HW 7 part 2:", {s: sp.simplify(v) for s, v in sol.items()})

# HW 7 part 3 (problem 10.25): flexibility matrices
print("HW 7 part 3:", sp.Matrix([[24, 3 * L], [3 * L, L**2]]).inv(),
      sp.Matrix([[24, 0, 3 * L], [0, 2 * L**2, L / 2], [3 * L, L / 2, L**2]]).inv())

# HW 9: shear and moment diagrams of the continuous beam
segments = [  # (x0, x1, V(x)) from the MATLAB piecewise definition
    (0, 1, lambda x: 7.62e4 - 100e3 * x), (1, 2, lambda x: -23800 + 0 * x), (2, 3, lambda x: -27950 + 0 * x),
    (3, 4, lambda x: 22050 + 0 * x), (4, 5, lambda x: 22050 - 25e3 * (x - 4)), (5, 6, lambda x: 11250 - 25e3 * (x - 5)),
    (6, 7, lambda x: -13750 + 0 * x),
]
x = np.concatenate([np.linspace(a, b, 100) for a, b, _ in segments])
V = np.concatenate([f(np.linspace(a, b, 100)) for a, b, f in segments])
from scipy.integrate import cumulative_trapezoid  # noqa: E402

M = cumulative_trapezoid(V, x, initial=0)
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.plot(x, V)
ax1.set(ylabel="Shear (N)", title="Shear diagram")
ax2.plot(x, M)
ax2.set(ylabel="Moment (N m)", xlabel="Distance (m)", title="Moment diagram")
for ax in (ax1, ax2):
    ax.grid(True)
# Reproduces the MATLAB numbers (its shear starts at 7.62e4 N although the reaction listed is 7.62e3 N).
x_max = 381 / 500  # dM/dx = 0 in the first span
M_max = -200 * x_max * (250 * x_max - 381)
sigma_max = M_max * 55e-3 / 2 / (36e-3 * 55e-3**3 / 12)
print(f"HW 9: maximum moment {M_max:.1f} N m at x = {x_max:.3f} m, maximum bending stress {sigma_max / 1e6:.3f} MPa")
plt.show()
