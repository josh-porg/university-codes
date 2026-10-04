"""DMD book chapter 3 (``Algorithm_3_1``): Koopman embedding of a nonlinear system with a slow manifold.

``x1' = mu x1``, ``x2' = lambda (x2 - x1^2)`` becomes linear in
``y = (x1, x2, x1^2)``: ``y' = [[mu, 0, 0], [0, lambda, -lambda], [0, 0, 2 mu]] y``.
Trajectories of the linear system are drawn with the attracting manifold
``y2 = y1^2``, the invariant set ``y3 = y1^2`` and the stable subspace of
the Koopman operator.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

mu, lam = -0.05, -1.0
A = np.array([[mu, 0, 0], [0, lam, -lam], [0, 0, 2 * mu]])
w, T = np.linalg.eig(A)
k = np.argmin(np.abs(w - 2 * mu))
slope = T[2, k] / T[1, k]  # slope of the stable subspace
print(f"Koopman eigenvalues {np.sort(w)}, stable-subspace slope y3/y2 = {slope:.4f}")

tspan = np.arange(0, 1000.01, 0.01)
trajs = [solve_ivp(lambda t, y: A @ y, (0, 1000), y0, t_eval=tspan, rtol=1e-8).y
         for y0 in ([1.5, -1, 2.25], [1, -1, 1], [2, -1, 4])]

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(projection="3d")
Xm, Zm = np.meshgrid(np.arange(-2, 2.01, 0.05), np.arange(-1, 4.01, 0.05))
ax.plot_surface(Xm, Xm**2, Zm, color="r", alpha=0.1)
Xm, Ym = np.meshgrid(np.arange(-2, 2.01, 0.05), np.arange(-1, 4.01, 0.05))
ax.plot_surface(Xm, Ym, Xm**2, color="b", alpha=0.1)
Xm, Ym = np.meshgrid(np.arange(-2, 2.01, 0.05), np.arange(0, 4.01, 0.05))
ax.plot_surface(Xm, Ym, slope * Ym, color=(0.3, 0.7, 0.3), alpha=0.7)
x = np.arange(-2, 2.01, 0.01)
ax.plot(x, x**2 / slope, x**2, "-g", lw=2)
ax.plot(x, x**2, x**2, "--r", lw=2)
ax.plot(x, x**2, -1 + 0 * x, "r--", lw=2)
for y in trajs:
    ax.plot(y[0], y[1], -1 + 0 * y[0], "k-", lw=1)
    ax.plot(y[0], y[1], y[2], "k", lw=1.5)
ax.plot([0, 0], [0, 0], [0, -1], "ko")
ax.set(xlim=(-4, 4), ylim=(-1, 4), zlim=(-1, 4), xlabel="y1", ylabel="y2", zlabel="y3")
ax.view_init(8, -15)
plt.show()
