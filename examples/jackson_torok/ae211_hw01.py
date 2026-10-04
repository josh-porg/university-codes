"""Jackson Torok, AE 211 homework 1 (``Torok_Jackson_HW1_MATLAB``): Moore, *MATLAB for Engineers*,
problems 2.3-2.17 and 3.4-3.23 (arithmetic, vectors, tables, statistics, complex impedance).

Fixes: problem 3.17 took ``size = length(g)`` of the gravity constant
instead of the grade vector ``G``; MATLAB's ``mode`` returns the smallest
of tied modes, kept here.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

rng = np.random.default_rng(0)

# 2.3
print("2.3:", 5**2, (5 + 3) / (5 * 6), np.sqrt(4 + 6**3), 9 + 6 / 12 + 7 * 5 ** (2 + 3),
      1 + 5 * 3 / 6**2 + 2 ** (2 - 4) / 1 / 5.5)
# 2.4
print(f"2.4: area(r=5) = {np.pi * 5**2:.4f}, surface(r=10) = {4 * np.pi * 10**2:.4f}, "
      f"volume(r=2) = {4 / 3 * np.pi * 2**3:.4f}")
# 2.8
h = np.array([1, 5, 12])
print("2.8: cylinder volumes", np.pi * 3**2 * h)
area_tri = 0.5 * np.array([2, 4, 6]) * 12
print("     triangle areas", area_tri, "prism volumes", area_tri * 6)
# 2.12
print("2.12:", np.linspace(4, 20, 15), np.logspace(1, 3, 10), sep="\n")
# 2.14
t = np.linspace(0, 100, 101)
print("2.14: time / distance (first rows)\n", np.column_stack([t, 0.5 * 9.8 * t**2])[:5])
# 2.17
print(f"2.17: Earth-Moon force = {6.673e-11 * 6e24 * 7.4e22 / 3.9e8**2:.4e} N")
# 3.4
T = np.linspace(100, 500, 9)
print("3.4: T, k\n", np.column_stack([T, 1200 / 60 * np.exp(-8000 / (1.987 * T))]))
# 3.10
theta = np.round(np.linspace(0, 2 * np.pi, 64), 1)
print("3.10: theta, sin, cos, tan (first rows)\n", np.column_stack([theta, np.sin(theta), np.cos(theta), np.tan(theta)])[:5])
# 3.17
G = np.array([68, 83, 61, 70, 75, 82, 57, 5, 76, 85, 62, 71, 96, 78, 76, 68, 72, 75, 83, 93])
print(f"3.17: n = {G.size}, sorted {np.sort(G)}, mean {G.mean():.2f}, median {np.median(G)}, "
      f"mode {stats.mode(G, keepdims=False).mode}, std {G.std(ddof=1):.4f}")
print("      the median is more representative: the outlier (5) pulls the mean down")
# 3.21
temps = rng.normal(70, 2, 241)
t = np.linspace(0, 120, 241)
plt.plot(t, temps)
print(f"3.21: max temp {temps.max():.3f} at {t[temps.argmax()]:g} min, min temp {temps.min():.3f} at {t[temps.argmin()]:g} min")
# 3.23
R, f, v, C, L = 5.0, 15e3, 10.0, 1e-9, 200e-3
w = 2 * np.pi * f
Z = 1 / (1j * w * C) + 1j * w * L + R
I = v / Z  # noqa: E741
print(f"3.23: |I| = {abs(I):.4e} A, angle = {np.angle(I):.4f} rad")
plt.show()
