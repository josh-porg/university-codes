"""AE 725 final exam problems 2, 4 and 5.

* Problem 2: min x1^2 + x2^2 inside a circle and left of x1 = 5 (``fmincon`` -> SLSQP).
* Problem 4: two-element cantilever, tip displacement and its sensitivity to
  I1 and I2 by the direct and adjoint methods (``Final_Exam_Problem_4_FEA``;
  the symbolic ``Final_Exam_Problem_4`` gives the same matrices).
* Problem 5: one CSD iteration - QP subproblem, descent function and Armijo step.
  The MATLAB typed the linearised constraint vector as b = [-5, -23]; with
  g2 = 2 - 25 = -23 it is b = -g = [-5, 23].
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

from unicodes.optimize import descent_function, inexact_step_size, qp_direction
from unicodes.structures import assemble_beam, displacement_sensitivity

# Problem 2. The MATLAB plotted the contours of x1^2 x2^2 but minimised x1^2 + x2^2.
cons = [{"type": "ineq", "fun": lambda x: 6 - (x[0] - 6) ** 2 - (x[1] - 6) ** 2},
        {"type": "ineq", "fun": lambda x: 5 - x[0]}]
r = minimize(lambda x: x @ x, [4.0, 4.0], method="SLSQP", constraints=cons)
print("Problem 2 optimum:", r.x, "f =", r.fun)
X, Y = np.meshgrid(np.linspace(-10, 10, 300), np.linspace(-10, 10, 300))
fig, ax = plt.subplots()
ax.contourf(X, Y, X**2 + Y**2, 15, cmap="summer_r")
ax.contour(X, Y, (X - 6) ** 2 + (Y - 6) ** 2 - 6, [0], colors="b")
ax.axvline(5, color="m")
ax.plot(*r.x, "*", color="purple", ms=12)
ax.set(xlabel="$x_1$", ylabel="$x_2$", title="Final exam problem 2")

# Problem 4: E = 4, L = 1, I1 = 2, I2 = 1, load P = 36 at the tip, clamped at node 1.
E, L, I1, I2, P = 4.0, 1.0, 2.0, 1.0, 36.0
K = assemble_beam(E, [I1, I2], L)
free = np.arange(2, 6)
F = np.zeros(6)
F[4] = P
u = np.linalg.solve(K[np.ix_(free, free)], F[free])
print("Problem 4 free displacements [w2, theta2, w3, theta3]:", u)
dK1 = assemble_beam(E, [1.0, 0.0], L)  # dK/dI1
dK2 = assemble_beam(E, [0.0, 1.0], L)
tip = np.array([0, 0, 1, 0])
for name, dK in (("I1", dK1), ("I2", dK2)):
    direct = displacement_sensitivity(K, dK, F, free, tip)
    adjoint = displacement_sensitivity(K, dK, F, free, tip, adjoint=True)
    print(f"  d w3 / d {name}: direct {direct:.4f}, adjoint {adjoint:.4f}")

# Problem 5: f = x1^2 + x2^2 - 2x1 - 2x2, g1 = x1 + x2 - 4, g2 = 2 - x1^2 at x = (5, 4)
x = np.array([5.0, 4.0])


def fun(x):
    return x[0] ** 2 + x[1] ** 2 - 2 * x[0] - 2 * x[1], np.array([2 * x[0] - 2, 2 * x[1] - 2])


def constraints(x):
    return np.array([x[0] + x[1] - 4, 2 - x[0] ** 2]), np.array([[1, 1], [-2 * x[0], 0]])


f0, c = fun(x)
g, A = constraints(x)
d, u = qp_direction(c, A, -g)
print("Problem 5: c =", c, "A =", A.tolist(), "b =", -g)
print("  QP direction d =", d, "multipliers u =", u)
phi, R = descent_function(fun, constraints, 1.0, u)
alpha = inexact_step_size(phi, x, d, 0.5, 0.5)
print(f"  R = {R:g}, alpha = {alpha:g}, x_new = {x + alpha * d}, f_new = {fun(x + alpha * d)[0]:.4f}")

t = np.linspace(0, 1, 200)
fig, ax = plt.subplots()
ax.plot(t, [phi(x + ti * d) for ti in t], "b", label=r"$\phi(x + t d)$")
ax.plot(t, phi(x) - t * 0.5 * (d @ d), "--r", label="Armijo line")
ax.axvline(alpha, color="k", label="accepted step")
ax.set(xlabel="t", title="Final exam problem 5 line search")
ax.legend()

plt.show()
