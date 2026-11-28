"""Jackson Torok, AE 211 exams 1 and 2 and the final (``Torok__Jackson_Exam1``,
``Torok_Jackson_Exam2``, ``exam``): the computational questions.

Written answers (good practice, variable naming rules, operator
precedence...) are not repeated here. Exam 1 Q10 noted that x + y + z of
the given system does not equal the 6 the question expected.
The final's radius prompt is ``--radius``.
"""

import argparse
import math

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad, solve_ivp, trapezoid
from scipy.special import erf

p = argparse.ArgumentParser()
p.add_argument("--radius", type=float, default=2.0)
a = p.parse_args()

print("=== Exam 1")
A, B = np.array([[2, 6], [5, 4]]), np.array([[7, 9], [1, 2]])
print("Q5: matrix product\n", A @ B, "\n element-wise\n", A * B)
aa = np.array([[1, 2, 3, 9, 8, 7], [4, 5, 6, 6, 5, 4], [7, 8, 9, 3, 2, 1], [4, 5, 6, 6, 5, 4], [7, 8, 9, 3, 2, 1],
               [1, 2, 3, 9, 8, 7]])
print("Q6:", aa[:5, 0], "\n", aa[:5, 2:6])
print(f"Q7: P(four of a kind in 4 cards) = {100 * 13 / math.comb(52, 4):.6f} %")
v1, v2 = np.array([3, 2, 4]), np.array([1, 3, 4])
print(f"Q8: angle = {np.degrees(np.arccos(v1 @ v2 / (np.linalg.norm(v1) * np.linalg.norm(v2)))):.4f} deg")
a9 = np.array([[2, 8, 15], [25, 9, 32], [14, 28, 51]])
b9, c9 = np.array([16, 17, 38]), np.array([7, 19, 34, 12])
d9 = a9[:, 2]
e9 = b9 @ d9
print(f"Q9: d = {d9}, e = {e9}, f = {np.cross(b9, d9)}, g = {np.sum(b9 + d9)}, i = {e9 + c9[2] + a9[1, 1]}, "
      f"det = {np.linalg.det(a9):.4f}")
sol = np.linalg.solve([[3, 4, -5], [-1, 2, 1], [2, -3, 4]], [-4, 6, 8])
print(f"Q10: x, y, z = {sol}, sum = {sol.sum():.4f}")
g, v = 9.8, 100.0
print(f"Q11: flight time at 40 deg = {2 * v * np.sin(np.radians(40)) / g:.4f} s")
thet = np.arange(91)
rng = v * np.sin(np.radians(thet)) / g * v * np.cos(np.radians(thet))
print(f"     the MATLAB's 'range' (half the true range) peaks at {rng.max():.2f} m at {thet[rng.argmax()]} deg")
plt.figure(num="Exam 1 Q11")
for ang in (25, 45, 75):
    t = np.linspace(0, 2 * v * np.sin(np.radians(ang)) / g, 200)
    plt.plot(t * v * np.cos(np.radians(ang)), t * v * np.sin(np.radians(ang)) - 0.5 * g * t**2, label=f"{ang} degrees")
plt.legend()

print("=== Exam 2")
X, Y = np.meshgrid(np.linspace(-2, 2, 100), np.linspace(-2, 2, 100))
Z = np.sin(3 * X + 2 * Y) * np.exp(-X**2 - Y**2)
fig = plt.figure(num="Exam 2 Q3")
fig.add_subplot(2, 2, 1, projection="3d").plot_wireframe(X, Y, Z, rstride=5, cstride=5)
fig.add_subplot(2, 2, 2, projection="3d").plot_surface(X, Y, Z)
ax = fig.add_subplot(2, 2, 3, projection="3d")
ax.plot_surface(X, Y, Z)
ax.contour(X, Y, Z, offset=-1)
fig.add_subplot(2, 2, 4).contour(X, Y, Z)
students = np.array([[2001, 20, 19], [2002, 23, 24], [2003, 19, 25], [2004, 30, 29], [2005, 25, 30], [2006, 15, 40],
                     [2007, 32, 29]])
for year, n_in, n_out in students[students[:, 1] > students[:, 2]]:
    print(f"Q4: in year {year}, there are {n_in} in-state students and {n_out} out-of-state students")
x = np.linspace(0, 10, 40)
f = np.select([x < 0, x <= 6, x <= 8], [0, x**2, 36], 4.5 * x)
plt.figure(num="Exam 2 Q5")
plt.plot(x, f)
guess, it = 0.0, 0
while abs(np.cos(np.pi / 6) - guess) > 1e-5:
    guess += (-1) ** it * (np.pi / 6) ** (2 * it) / math.factorial(2 * it)
    it += 1
print(f"Q6: cos(pi/6) series = {guess:.6f}, error {np.cos(np.pi / 6) - guess:.3e}, {it} iterations")

print("=== Final exam")
print(f"Q3: area = {np.pi * a.radius**2:.4f}")
a6, b6, x6 = -3, 2, 0
i6 = b6 / 2
print("Q6:", (a6 <= x6 and b6 >= 2 and i6 > 1) or (b6 * x6 != 0))
A7 = np.array([[2, 4, 6], [1, 3, 5], [7, 8, 9]])
C7 = np.arange(1, 4)
print(f"Q7: D{{1}}(2,3) + D{{3}}(2) = {A7[1, 2] + C7[1]}")
x8 = 0
print("Q8: f(0) =", -x8**2 if x8 < 0 else x8**2 if x8 <= 2 else x8**3 if x8 <= 6 else x8)
print("Q9:\n", np.pi * np.ones((4, 6)))
ode = solve_ivp(lambda t, y: t**2 + y, (0, 1), [0.0], rtol=1e-8)
print(f"Q10: y' = t^2 + y, y(0) = 0: y(1) = {ode.y[0, -1]:.6f} (exact {2 * np.e - 5:.6f})")


def fibonacci(first, second, num):
    f = [first, second]
    while len(f) < num:
        f.append(f[-1] + f[-2])
    return np.array(f)


a12, b12 = fibonacci(1, 1, 10), fibonacci(5, 8, 10)
print("Q12:", a12, b12, a12 * b12)
q = 2
while True:
    xs = np.linspace(-1, 1, q)
    pi_guess = 2 * trapezoid(np.sqrt(1 - xs**2), xs)
    if np.pi - pi_guess <= 1e-6:
        break
    q += 1
print(f"Q13: trapezoidal pi = {pi_guess:.8f} with {q} points")
x = np.linspace(0, 3, 100)
errf = np.array([2 / np.sqrt(np.pi) * quad(lambda s: np.exp(-s**2), 0, xi)[0] for xi in x])
print(f"Q14: max |quad erf - scipy erf| = {np.abs(errf - erf(x)).max():.2e}")
fig, ax = plt.subplots(1, 2, num="Final Q14")
ax[0].plot(x, errf, x, np.gradient(errf, x))
ax[0].set(title="Error Function", xlabel="Values of x", ylabel="erf(x)")
ax[1].plot(x, errf, x, np.polyval(np.polyfit(x, errf, 3), x))
ax[1].legend(["erf(x)", "Polynomial Approx"], loc="lower right")
plt.show()
