"""Jackson Torok, AE 211 homework 5 (``Torok_Jackson_HW5_MATLAB`` with ``num_grain``, ``distance``,
``height``, ``plotndfhs``): user-defined functions (6.1-6.14).

ASTM grain count ``2^(n - 1)``, horizon distance ``sqrt(2 r h + h^2)`` on
Earth and Mars, and a rocket height ``h(t) = -4.9 t^2 + 125 t + 500``
(NaN for t < 0) with its apex and ground-impact times. ``plotndfhs`` (a
parabola test plot) is the last figure.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq


def num_grain(n):
    return 2.0 ** (n - 1)


def distance(r, h):
    """Horizon distance table: rows are heights, columns radii (the MATLAB meshgrid)."""
    r, h = np.meshgrid(r, h)
    return np.sqrt(2 * r * h + h**2)


def height(t):
    t = np.asarray(t, float)
    return np.where(t >= 0, -9.8 / 2 * t**2 + 125 * t + 500, np.nan)


# 6.1
n = np.arange(10, 101)
plt.semilogy(n, num_grain(n))
plt.title("Grains per square inch at x100 vs ASTM grain size")
# 6.6
h = np.arange(0, 10001, 100) / 5280
dist = distance([7926, 4217], h)
print("6.6: horizon distance (mi) on Earth / Mars, every 1000 ft:\n", np.column_stack([h[::10] * 5280, dist[::10]]))
# 6.7
t = np.arange(0, 30.01, 0.5)
plt.figure()
plt.plot(t, height(t))
print(f"6.7: maximum height at t = {t[np.nanargmax(height(t))]:.1f} s")
# 6.14
plt.figure()
tt = np.linspace(0, 60, 400)
plt.plot(tt, height(tt))
print(f"6.14: hits the ground at t = {brentq(height, 1, 60):.4f} s")
# plotndfhs
plt.figure()
x = np.arange(-10, 10.01, 0.1)
plt.plot(x, x**2)
plt.show()
