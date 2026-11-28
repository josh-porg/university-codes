"""AE 746 project 1: periodic linear advection with three schemes and a grid-convergence study
(``Project_1``, ``Project_1_part_2``).

``u_t + u_x = 0`` on [0, 10], CFL 0.25, initial condition ``1 - cos x`` on
[0, 2 pi] (and 1 on [7, 9] in part 1), solution at t = 20.

Fixes: part 1 used ``dt = sigma / (c dx)`` (part 2 corrected it to
``sigma dx / c``); the square pulse test read ``x <= 7 && x <= 9``.
"""

import matplotlib.pyplot as plt
import numpy as np

from unicodes.cfd.advection import SCHEMES, advect

L, c, sigma, t_end = 10.0, 1.0, 0.25, 20.0


def initial(x, pulse=True):
    u = np.where(x <= 2 * np.pi, 1 - np.cos(x), 0.0)
    return np.where(pulse & (x >= 7) & (x <= 9), 1.0, u)


# Part 1: compare the schemes
for n in (40, 80):
    dx = L / n
    x = dx * np.arange(n)
    steps = int(round(t_end * c / (sigma * dx)))
    fig, ax = plt.subplots()
    for name in SCHEMES:
        u = advect(initial(x), sigma, steps, name)
        ax.plot(x, np.clip(u, -3, 5), label=name)
    ax.plot(x, initial((x - c * t_end) % L), "k--", label="exact")
    ax.set(title=f"Solution at t = 20, dx = {dx}", xlabel="x", ylabel="u")
    ax.legend(loc="lower left")

# Part 2: convergence of second-order upwind / SSP-RK2 (smooth initial condition only)
errors = []
fig, ax = plt.subplots()
for n in (40, 80, 160):
    dx = L / n
    x = dx * np.arange(1, n + 1)
    steps = int(round(t_end / (sigma * dx / c)))
    u = advect(initial(x, pulse=False), sigma, steps)
    exact = initial((x - t_end) % L, pulse=False)
    errors.append(np.sqrt(np.mean((u - exact) ** 2)))
    ax.plot(x, u, label=f"{n} grid points")
ax.plot(x, exact, "k--", label="exact")
ax.set(title="O(2) upwind SSP-RK2 at t = 20", xlabel="x", ylabel="u")
ax.legend()
errors = np.array(errors)
# The MATLAB used the sum of squared errors, whose log2 ratio is 2p - 1, not p.
print("RMS errors:", errors, "observed orders:", np.log2(errors[:-1] / errors[1:]))
plt.show()
