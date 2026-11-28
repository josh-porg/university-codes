"""AE 507 aerospace structures: homework 1, 2.2.4 and 2.5, the 2/9/2022 in-class problem and the
design-project wall thickness (``HW_1``, ``HW_2_2_4``, ``HW_2_5``, ``in_class_2_9_2022``,
``Design_project_thickness_solver``).
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq

# HW 1: centripetal and tangential acceleration
w = np.array([0, -0.6857142857, 0])
r = np.array([-20.83333333333, 0, -2.032520324])
alpha = np.array([0, 13.14, 0])
wwr = np.cross(w, np.cross(w, r))
print("HW 1: w x (w x r) =", wwr, "|.| =", np.linalg.norm(wwr), " alpha x r =", np.cross(alpha, r))

# HW 2.2.4: reactions of a strut-braced fitting (moment and force balance)
N = 2 * np.array([3000.0, 12500, -1000])
d1, d2, d3, d4, d5 = 12, 24, 8, 26, 30
rN, rB, rC, rD = np.array([0, -d5 - d3, 0]), np.array([0, d2, 0]), np.array([d4, d2, -d1]), np.array([d4, d2, d1])
eC, eD = rC / np.linalg.norm(rC), rD / np.linalg.norm(rD)
# unknowns: Bx, By, Bz, c, d  -> sum F = 0 (3) and sum M about the origin = 0 (3, one redundant)
A = np.zeros((6, 5))
A[:3, :3] = np.eye(3)
A[:3, 3], A[:3, 4] = eC, eD
skew = lambda v: np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])  # noqa: E731
A[3:, :3] = skew(rB)
A[3:, 3], A[3:, 4] = np.cross(rC, eC), np.cross(rD, eD)
b = -np.r_[N, np.cross(rN, N)]
x = np.linalg.lstsq(A, b, rcond=None)[0]
print("HW 2.2.4: B =", np.round(x[:3], 2), " |C| =", round(x[3], 2), " |D| =", round(x[4], 2),
      " residual", np.linalg.norm(A @ x - b))

# HW 2.5: shear, moment, slope and deflection of a wing with uniform then tapering lift
c1, V_root, M_root = 8.42, 1137.15, 76867.9  # lbf/in, lbf, lbf in
xs = np.linspace(15, 160, 146)
lift = np.where(xs <= 140, c1, c1 + (0 - c1) / 20 * (xs - 140))
V = V_root - cumulative_trapezoid(lift, xs, initial=0)
M = M_root + cumulative_trapezoid(V, xs, initial=0)  # sign convention of the MATLAB
EI = 1.0  # the MATLAB left E I = 1 as a TODO; deflections below are per unit E I
theta = cumulative_trapezoid(M / EI, xs, initial=0)
v = cumulative_trapezoid(theta, xs, initial=0)
fig, axes = plt.subplots(5, 1, sharex=True, figsize=(6, 9))
for ax, y, lab in zip(axes, (lift, V, M, theta, v), ("lift (lbf/in)", "V (lbf)", "M (lbf in)", "EI theta", "EI v")):
    ax.plot(xs, y)
    ax.set_ylabel(lab)
axes[-1].set_xlabel("x (in)")
print(f"HW 2.5: tip shear {V[-1]:.2f} lbf, tip moment {M[-1]:.1f} lbf in")

# In class 2/9/2022: wall thickness of a 1.5 in square tube with I = 600 * 0.75 / 7000
I_req = 600 * 0.75 / 7000
t = brentq(lambda t: (1.5**4 - (1.5 - 2 * t) ** 4) / 12 - I_req, 1e-6, 0.75)
print(f"In class: t = {t:.5f} in")

# Design project: square-tube wall thickness for a tip-rotation margin of safety of 0.1.
# The MATLAB wrote the second moment as (s^4 - (s - 2t_v)(s^4 - (s - 2t_h)^3)) / 12; the
# square tube I = (s^4 - (s - 2t)^4) / 12 with t_h = t_v is used here.
theta_all, E, s, L, W, n, F_lm, F_gm, MS = np.deg2rad(5), 4e6, 2.0, 18.0, 0.277, 3.431, 24.0, 1.0, 0.1


def margin(t):
    I = (s**4 - (s - 2 * t) ** 4) / 12
    rotation = W / L * n * L**3 / (6 * E * I) + (F_lm - F_gm * n) * L**2 / (2 * E * I)
    return theta_all / rotation - 1 - MS


t = brentq(margin, 1e-5, s / 2)
print(f"Design project: wall thickness {t:.4f} in")
plt.show()
