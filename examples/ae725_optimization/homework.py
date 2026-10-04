"""AE 725 homework: graphical optimisation, KKT checks, fmincon problems and CSD.

* HW 1 problem 3.51 - water-tower column (``Homework_1_Problem_3_51``)
* HW 3 problems 4.43 and 4.59 - KKT cases by linear solves
* HW 5 problems 5.4, 5.62, 7.5, 7.10 (``fmincon`` -> SLSQP)
* HW 7 problems 10.52 (steepest descent with exact line search) and 10.57 (Hessian)
* HW 10 problem 13.1 - beam design by constrained steepest descent

The constraint regions are shaded from the constraint functions themselves;
the MATLAB drew them from hand-solved boundary curves (several of which
were only approximate). ``Homework_5_Problem_SP1`` did not run (undefined
variables, a syntax error) and is not ported; the symbolic multiplier
derivations are replaced by numerical KKT checks.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

from unicodes.optimize import constrained_steepest_descent, golden_section_step


def shade(ax, X, Y, g, color, label):
    ax.contourf(X, Y, (g > 0).astype(float), levels=[0.5, 1.5], colors=[color], alpha=0.3)
    ax.contour(X, Y, g, levels=[0], colors=[color])
    ax.plot([], [], color=color, label=label)


def slsqp(f, x0, g, bounds=None):
    """``fmincon`` with ``g(x) <= 0``."""
    return minimize(f, x0, method="SLSQP", bounds=bounds, constraints={"type": "ineq", "fun": lambda x: -np.asarray(g(x))},
                    options={"ftol": 1e-12, "maxiter": 500})


# HW 1 problem 3.51 / HW 5 problem 7.5: water tower column, x = [d_o, t]
h, H, D, w, e = 10.0, 30.0, 10.0, 700.0, 0.1
gamma_w, gamma_s, E, t_t = 1e4, 80e3, 210e9, 1.5e-2
Delta, sigma_b = 2e-2, 165e6
V = 1.2 * np.pi * D**2 * h
P = V * gamma_w + t_t * 1.25 * np.pi * D**2 * gamma_s
W = w * 2 / 3 * D * h
c_delta = W * H**2 / (12 * E) * (4 * H + 3 * h) + H / (2 * E) * (0.5 * W * h + P * e) * (H + h)


def tower_constraints(x):
    d_o, t = x
    A = np.pi * t * (d_o - t)
    I = np.pi / 64 * (d_o**4 - (d_o - 2 * t) ** 4)
    r = np.sqrt(I / A)
    sigma_a = 12 * np.pi**2 * E / (92 * H / r)
    delta = c_delta / I
    M = W * (H + 0.5 * h) + (delta + e) * P
    return np.array([
        (d_o - t) / 2 - 2, 0.35 - (d_o - t) / 2, t - 0.2, 0.01 - t, d_o / t - 92,
        delta - Delta,
        # The MATLAB had f_a/sigma_a - f_b/sigma_b; the interaction check is f_a/sigma_a + f_b/sigma_b - 1.
        P / A / sigma_a + M / (2 * I) * d_o / sigma_b - 1,
    ])  # fmt: skip


res = slsqp(lambda x: x[1] * (x[0] - x[1]), [1.0, 0.05], tower_constraints, bounds=[(0.1, 99), (1e-3, 99)])
print(f"HW 7.5 water tower: d_o = {res.x[0]:.4f} m, t = {res.x[1]:.4f} m, f = {res.fun:.5f}")
d_o, t = np.meshgrid(np.linspace(0.05, 5, 300), np.linspace(1e-3, 0.3, 300))
with np.errstate(all="ignore"):
    G = np.array([[tower_constraints((a, b)) for a, b in zip(ra, rb)] for ra, rb in zip(d_o, t)])
fig, ax = plt.subplots()
ax.contourf(d_o, t, t * (d_o + t), 10, cmap="winter")
for i, col in enumerate(plt.cm.autumn(np.linspace(0, 1, G.shape[-1]))):
    shade(ax, d_o, t, G[..., i], col, f"g{i + 1}")
ax.plot(*res.x, "y*", ms=12)
ax.set(xlabel="outer diameter $d_o$ (m)", ylabel="thickness $t$ (m)", title="HW 1 3.51 / HW 5 7.5")
ax.legend(fontsize=7)

# HW 3 problem 4.43: min 4x1^2 + 3x2^2 + 5x1x2 - 8x1 s.t. x1 + x2 = 4
K = np.array([[8, 5, 1], [5, 6, 1], [1, 1, 0]])
print("HW 4.43 [x1, x2, v] =", np.linalg.solve(K, [8, 0, 4]))

# HW 3 problem 4.59: KKT switching cases, unknowns [x1, x2, u1, u2]
A = np.array([[2, 0, -1, -1], [0, 2, -1, 1], [-1, -1, 0, 0], [-1, 1, 0, 0]], float)
c = np.array([2, 2, -4, -2], float)
for name, idx in [("s1 = s2 = 0", [0, 1, 2, 3]), ("u1 = s2 = 0", [0, 1, 3]), ("s1 = u2 = 0", [0, 1, 2]),
                  ("u1 = u2 = 0", [0, 1])]:
    print(f"HW 4.59 {name}:", np.linalg.solve(A[np.ix_(idx, idx)], c[idx]))

# HW 5 problem 5.4: min 4x1^2 + 3x2^2 - 5x1x2 - 8x1 s.t. x1 + x2 = 4.
# The KKT solution is (13/6, 11/6); the MATLAB marked (41/19, 35/19) from a mistyped KKT matrix.
f54 = lambda x: 4 * x[0] ** 2 + 3 * x[1] ** 2 - 5 * x[0] * x[1] - 8 * x[0]  # noqa: E731
r = minimize(f54, [0, 0], method="SLSQP", constraints={"type": "eq", "fun": lambda x: x[0] + x[1] - 4})
print("HW 5.4:", r.x, "expected", [13 / 6, 11 / 6])

# HW 5 problem 5.62: thin-walled column, x = [R, t]. Buckling and R/t = 50 active give
# R = (200 P l^2 / (pi^3 E))^(1/4); the MATLAB used 100.
P62, l, E62, sa = 50e3, 5.0, 210e9, 250e6
beta = (200 * P62 * l**2 / (np.pi**3 * E62)) ** 0.25
g62 = lambda x: [P62 / (2 * np.pi * x[0] * x[1]) / sa - 1, 1 - np.pi**3 * E62 * x[0] ** 3 * x[1] / (4 * l**2) / P62,  # noqa: E731
                 x[0] / x[1] / 50 - 1]
r = slsqp(lambda x: x[0] * x[1] * 1e4, [0.1, 0.005], g62, bounds=[(1e-3, 1), (1e-5, 0.01)])
print(f"HW 5.62: R, t = {r.x}, analytic ({beta:.4f}, {beta / 50:.5f})")

# HW 5 problem 7.10: channel section, x = [t1, t2, b, h]
P710, E710, theta, L = -70e3, 200e9, np.deg2rad(45), 1.5


def channel_constraints(x):
    t1, t2, b, h = x
    A = b * (2 * t1 - t2)
    Q = b * t1 * (t1 + h) / 2 + t2 * h**2 / 8
    Iy = (b * (h + 2 * t1) ** 3 - (b - t2) * h**3) / 12
    Iz = (2 * t1 * b**3 + h * t2**3) / 12
    Px, v = P710 * np.cos(theta), P710 * np.sin(theta)
    sigma = Px * L * (h / 2 + t1) / Iy + Px / A
    tau = v * Q / (Iy * t2)
    delta = v * L**3 / (3 * E710 * Iy)
    # Magnitudes are used: with the compressive (negative) load the MATLAB constraints were always satisfied.
    return np.array([abs(sigma) / 100e6 - 1, abs(tau) / 60e6 - 1, abs(delta) / 0.015 - 1,
                     abs(P710) / (np.pi * E710 * Iy / (4 * L**2)) - 1, abs(P710) / (np.pi * E710 * Iz / (4 * L**2)) - 1,
                     0.1 - b, t1 - 0.01, t2 - 0.015, h - 0.15])


r = slsqp(lambda x: 2 * x[2] * x[0] + x[3] * x[1], [0.01, 0.015, 0.2, 0.15], channel_constraints,
          bounds=[(1e-4, 0.02), (1e-4, 0.02), (0, 99), (1e-3, 0.2)])
print("HW 7.10 [t1, t2, b, h] =", r.x, "area", r.fun, "max g", channel_constraints(r.x).max())

# HW 7 problem 10.52: two steepest-descent iterations with exact line search
f1052 = lambda x: x[0] ** 2 + 2 * x[1] ** 2 - 4 * x[0] - 2 * x[0] * x[1]  # noqa: E731
grad1052 = lambda x: np.array([2 * x[0] - 4 - 2 * x[1], 4 * x[1] - 2 * x[0]])  # noqa: E731
x = np.array([1.0, 1.0])
for i in range(2):
    d = -grad1052(x)
    alpha = golden_section_step(f1052, x, d, tol=1e-10)
    x = x + alpha * d
    print(f"HW 10.52 iteration {i + 1}: alpha = {alpha:.6f}, x = {x}, f = {f1052(x):.6f}")

# HW 7 problem 10.57: Hessian of x1^2 + 2x2^2 + 2x3^2 + 2x1x2 + 2x2x3
Hs = np.array([[2, 2, 0], [2, 4, 2], [0, 2, 4]])
print("HW 10.57 eigenvalues:", np.linalg.eigvalsh(Hs), "H @ [4, 8, 8] =", Hs @ [4, 8, 8])

# HW 10 problem 13.1: min b d (beam), x = [b, d] with 10 <= b, d <= 1000


def beam_fun(x):
    return x[0] * x[1], np.array([x[1], x[0]])


def beam_constraints(x):
    b, d = x
    g = np.array([2.4e8 / (b * d**2) - 10, 4.5e5 / (2 * b * d) - 2, d - 2 * b,
                  10 - b, 10 - d, b - 1e3, d - 1e3])  # bounds as constraints (quadprog lb/ub in the MATLAB)
    A = np.array([[-2.4e8 / (b**2 * d**2), -4.8e8 / (b * d**3)],
                  [-4.5e5 / (2 * b**2 * d), -4.5e5 / (2 * b * d**2)],
                  [-2, 1], [-1, 0], [0, -1], [1, 0], [0, 1]])
    # The MATLAB gradient of g1 w.r.t. d had an extra factor 1/2.
    return g, A


r = constrained_steepest_descent(beam_fun, beam_constraints, [250, 300], max_iterations=10000, R=1)
print(f"HW 13.1: b = {r.x[0]:.2f} mm, d = {r.x[1]:.2f} mm, f = {r.fun:.1f} after {r.iterations} iterations")

plt.show()
