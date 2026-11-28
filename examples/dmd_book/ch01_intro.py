"""DMD book chapter 1 (``Algorithm_1_2`` - ``Algorithm_1_5``, ``DMDfull``): DMD of two mixed
spatio-temporal patterns, compared with PCA and ICA.

``f1 = sech(x + 3) e^{2.3 i t}`` and ``f2 = sech(x) tanh(x) 2 e^{2.8 i t}``
on x in [-10, 10], t in [0, 4 pi]. Rank-2 DMD separates them with the
right frequencies; the first two PCA modes and ICA components (of the real
part) mix them. ``DMDfull`` (DMD with optional shift-stacking) is
:class:`unicodes.decomposition.DMD` with
:func:`~unicodes.decomposition.time_delay_stack`.
"""

import matplotlib.pyplot as plt
import numpy as np

from unicodes.decomposition import DMD, fast_ica, time_delay_stack

xi = np.linspace(-10, 10, 400)
t = np.linspace(0, 4 * np.pi, 200)
dt = t[1] - t[0]
Xg, T = np.meshgrid(xi, t)
f1 = 1 / np.cosh(Xg + 3) * np.exp(1j * 2.3 * T)
f2 = 1 / np.cosh(Xg) * np.tanh(Xg) * 2 * np.exp(1j * 2.8 * T)
f = f1 + f2
X = f.T

# Algorithm 1.3: rank-2 DMD
dmd = DMD(X, rank=2, dt=dt, truncation=None)
print("DMD omega:", np.round(dmd.omega, 6), "(true 2.3i and 2.8i)")
X_dmd = dmd.reconstruct()
print(f"relative reconstruction error {np.linalg.norm(X_dmd - X) / np.linalg.norm(X):.2e}")

# Algorithm 1.4: PCA
U, S, Vh = np.linalg.svd(X)
pc1, pc2 = U[:, 0], U[:, 1]

# Algorithm 1.5: ICA of the real part
ic, ict = fast_ica(X.real.T, 2, rng=0)  # spatial components (2, 400), time courses (200, 2)

# DMDfull with two stacks gives the same two frequencies
stacked = DMD(time_delay_stack(X, 2), rank=2, dt=dt, truncation=None)
print("DMD with 2 stacks:", np.round(stacked.omega, 6))

fig = plt.figure(figsize=(10, 7))
for k, (data, title) in enumerate(((f1, "f1"), (f2, "f2"), (f, "f = f1 + f2"), (X_dmd.T, "rank-2 DMD"))):
    ax = fig.add_subplot(2, 2, k + 1, projection="3d")
    ax.plot_surface(Xg, T, data.real, cmap="gray", linewidth=0)
    ax.view_init(60, -20)
    ax.set_title(title)
fig, ax = plt.subplots(3, 1, figsize=(7, 7), sharex=True)
ax[2].set_xlabel("x")
ax[0].plot(xi, np.abs(dmd.modes[:, 0]), xi, np.abs(dmd.modes[:, 1]))
ax[0].set_title("|DMD modes|")
ax[1].plot(xi, pc1.real, xi, pc2.real)
ax[1].set_title("PCA modes")
ax[2].plot(xi, ic[0], xi, ic[1])
ax[2].set_title("ICA modes (real part)")
plt.show()
