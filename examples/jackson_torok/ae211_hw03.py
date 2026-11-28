"""Jackson Torok, AE 211 homework 3 (``Torok_Jackson_HW3_MATLAB``): vector angles, linear systems,
centre of mass, moments, determinants and inverses (10.14-10.18), function plots (5.1, 5.4).

Fix: problem 10.14 printed ``inv(A)`` for all three matrices; here each
invertible matrix is inverted (B is singular).
"""

import matplotlib.pyplot as plt
import numpy as np

# Question 1
F1, F2 = np.array([10, 5, 6]), np.array([-8, 3, 9])
cos_t = np.clip(F1 @ F2 / (np.linalg.norm(F1) * np.linalg.norm(F2)), -1, 1)
print(f"Q1: angle between the forces = {np.degrees(np.arccos(cos_t)):.4f} deg")
# Question 2
sol = np.linalg.solve([[3, 4, 8], [-1, 1, 4], [-5, 4, 7]], [15, 23, 5])
print(f"Q2: x, y, z = {sol}, 7x - 2y + z = {7 * sol[0] - 2 * sol[1] + sol[2]:.4f} "
      "(differs from the expected value, as noted in the homework)")
# Question 3
m = np.array([15, 25, 10, 50, 30])
loc = np.array([[10, 15, 5], [15, 5, -10], [20, 25, 15], [5, 5, 5], [25, 35, 20]])
print("Q3: centre of mass", m @ loc / m.sum())
# Question 4
L, F, th = 10 / 12, np.array([750, 150, 500]), np.array([30, 120, 190])
print(f"Q4: moment about A = {L * np.sum(F * np.sin(np.radians(th - th[0]))):.4f} ft lb")
# 10.14
for name, M in (("A", [[2, -1], [2, 5]]), ("B", [[4, 2], [2, 1]]), ("C", [[2, 0, 0], [1, 2, 2], [5, -4, 0]])):
    det = np.linalg.det(M)
    print(f"10.14: det {name} = {det:.4f}" + (f", inverse =\n{np.linalg.inv(M)}" if abs(det) > 1e-12 else " (singular)"))
# 10.16
l, h, Fb, theta_f = 10 / 12, 5 / 12, 35, 55
theta_pos = 180 - np.degrees(np.arctan(h / l))
print(f"10.16: bracket moment = {np.hypot(l, h) * Fb * np.sin(np.radians(theta_pos - theta_f)):.4f} ft lb")
# 10.18
systems = (([[-2, 1], [1, 1]], [3, 10]), ([[5, 3, -1], [3, 2, 1], [4, -1, 3]], [10, 4, 12]),
           ([[3, 1, 1, 1], [1, -3, 7, 1], [2, 2, -3, 4], [1, 1, 1, 1]], [24, 12, 17, 0]))
for k, (A, b) in enumerate(systems):
    print(f"10.18{'abc'[k]}: solve {np.linalg.solve(A, b)}, inverse {np.linalg.inv(A) @ b}")
# 5.1
x = np.linspace(0, 10, 400)
fig, axes = plt.subplots(2, 2)
for ax, (f, title, lim) in zip(axes.flat, ((np.exp, "exp(x)", (0, 1e4)), (np.sin, "sin(x)", (-1, 1)),
                                            (lambda x: 5 * x**2 + 2 * x + 4, "5x^2+2x+4", (0, 500)),
                                            (np.sqrt, "sqrt(x)", (0, 4)))):
    ax.plot(x, f(x))
    ax.set(title=title, xlabel="X", ylabel="Y", xlim=(0, 10), ylim=lim)
    ax.grid(True)
# 5.4
x = np.linspace(-np.pi, np.pi, 400)
plt.figure()
plt.plot(x, np.sin(x), "--r", x, np.sin(2 * x), "-b", x, np.sin(3 * x), ":g")
plt.axis([-np.pi, np.pi, -1, 1])
plt.show()
