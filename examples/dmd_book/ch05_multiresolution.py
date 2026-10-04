"""DMD book chapter 5 (``mrDMD``, ``mrDMD_map``, ``mrDMD_demo``): multiresolution DMD of a noisy movie.

An 80 x 80 movie (dt = 0.01 s, 10 s) of three modes that switch on and
off: a Gaussian blob at 5.55 Hz (0-5 s), a square at 0.9 Hz (3-7 s) and a
large rectangle at 0.15 Hz (always), plus Gaussian noise (sigma 0.1).
Six levels of rank-10 mrDMD; the time-frequency map shows when each mode
is active. ``--animate`` plays the movie.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from unicodes.decomposition import mrdmd

p = argparse.ArgumentParser()
p.add_argument("--animate", action="store_true")
p.add_argument("--seed", type=int, default=0)
a = p.parse_args()
rng = np.random.default_rng(a.seed)

nx = ny = 80
T, dt = 10.0, 0.01
t = np.arange(1, int(round(T / dt)) + 1) * dt
Xg, Yg = np.meshgrid(np.arange(1, nx + 1), np.arange(1, ny + 1))
u2 = np.zeros((nx, ny))
u2[nx - 41:nx - 10, ny - 41:ny - 10] = 1
u3 = np.zeros((nx, ny))
u3[:nx - 20, :ny - 20] = 1
modes = [(np.exp(-((Xg - 40) ** 2 / 250 + (Yg - 40) ** 2 / 250)), 5.55, 1.0, (0, 5)),
         (u2, 0.9, 1.0, (3, 7)), (u3, 0.15, 0.5, (0, T))]
Xclean = np.zeros((nx * ny, t.size))
for ti in range(1, t.size + 1):
    snap = np.zeros((nx, ny))
    for u, f, A, (t0, t1) in modes:
        if round(t0 / dt) < ti < round(t1 / dt):
            snap += A * u * np.cos(2 * np.pi * f * dt * ti)
    Xclean[:, ti - 1] = snap.ravel(order="F")
X = Xclean + 0.1 * rng.standard_normal(Xclean.shape)

if a.animate:
    fig, ax = plt.subplots()
    im = ax.imshow(X[:, 0].reshape(nx, ny, order="F"), vmin=-np.ptp(X) / 2, vmax=np.ptp(X) / 2)
    for ti in range(0, t.size, 5):
        im.set_data(X[:, ti].reshape(nx, ny, order="F"))
        plt.pause(dt)

L, r = 6, 10
tree = mrdmd(X, dt, r, max_cycles=2, levels=L)
amp, low_f = tree.amplitude_map()
J = amp.shape[1]
for k in range(L):
    freqs = np.unique(np.round(np.concatenate([np.abs(w.omega.imag) for w in tree.level(k)]), 2))
    print(f"level {k}: cutoff {tree.level(k)[0].rho:.3f} Hz, slow-mode frequencies {freqs[freqs > 0][:8]}")
fig, ax = plt.subplots()
ax.imshow(-amp, aspect="auto", origin="lower", cmap="pink", extent=(0, T, 0, L))
ax.set_yticks(np.arange(L) + 0.5, [f"{np.floor(f * 10) / 10:g}" for f in low_f[:L]])
ax.set(xlabel="Time (sec)", ylabel="Freq. (Hz)")
ax.grid(True)
plt.figure()
plt.imshow(np.abs(tree.level(0)[0].modes[:nx * ny, 0]).reshape(nx, ny, order="F"))
plt.title("|slowest mode|, level 1")
plt.show()
