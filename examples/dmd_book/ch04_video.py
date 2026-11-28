"""DMD book chapter 4 (``Algorithm_4_1`` - ``Algorithm_4_4``): background/foreground separation.

A static background ``0.5 cos(x)`` plus a moving ``2 sech(x) tanh(x) e^{2.8 i t}``
foreground. DMD modes with ``|omega| < 1e-2`` are the background; each part
is reconstructed with its own least-squares amplitudes, as in the book.
The book asked for 50 modes of rank-3 data; singular values at round-off
level are dropped here (the 0.5 % truncation), otherwise the spurious modes
spoil the foreground amplitudes.
"""

import matplotlib.pyplot as plt
import numpy as np

from unicodes.decomposition import DMD

n, m = 200, 80
x = np.linspace(-15, 15, n)
t = np.linspace(0, 8 * np.pi, m)
dt = t[1] - t[0]
Xg, T = np.meshgrid(x, t)
f1 = 0.5 * np.cos(Xg) * (1 + 0 * T)
f2 = 1 / np.cosh(Xg) * np.tanh(Xg) * 2 * np.exp(1j * 2.8 * T)
X = (f1 + f2).T

dmd = DMD(X, rank=50, dt=dt)  # the data have rank 3: the truncation drops the round-off modes
bg = np.flatnonzero(np.abs(dmd.omega) < 1e-2)
fg = np.setdiff1d(np.arange(dmd.rank), bg)
print(f"{dmd.rank} modes kept; background modes {bg}, omega_bg = {np.round(dmd.omega[bg], 6)}")


def part(idx):
    b = np.linalg.lstsq(dmd.modes[:, idx], X[:, 0], rcond=None)[0]
    return dmd.modes[:, idx] @ (b[:, None] * np.exp(np.outer(dmd.omega[idx], t)))


X_bg, X_fg = part(bg), part(fg)
print(f"background error {np.linalg.norm(X_bg - f1.T) / np.linalg.norm(f1):.2e}, "
      f"foreground error {np.linalg.norm(X_fg - f2.T) / np.linalg.norm(f2):.2e}")
plt.figure()
plt.plot(dmd.omega.real, dmd.omega.imag, ".")
plt.title("DMD omega")
fig = plt.figure(figsize=(12, 4))
for k, (data, title) in enumerate(((X, "data"), (X_bg, "background"), (X_fg, "foreground"))):
    ax = fig.add_subplot(1, 3, k + 1, projection="3d")
    ax.plot_surface(Xg, T, data.T.real, cmap="gray", linewidth=0)
    ax.view_init(60, -20)
    ax.set_title(title)
plt.show()
