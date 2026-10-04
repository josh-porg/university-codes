"""DMD of two superposed spatio-temporal modes (``DMD_simple_functions``,
``DMD_benchmarking_againgst_simple_functions``, ``timeandSpacialWave``,
``timedynamics_experimentation``).

``f = sech(x + 3 - t) e^(2.3 i t) + 2 sech(x) tanh(x) e^(2.8 i t)`` (the
travelling version of the benchmarking script; ``--standing`` drops the
``-t``) on x in [-10, 10], t in [0, 4 pi]. The library DMD is checked
against a hand-written rank-2 exact DMD, as the MATLAB did for its
``DynamicModeDecomposer`` class.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from unicodes.decomposition import DMD

p = argparse.ArgumentParser()
p.add_argument("--standing", action="store_true")
p.add_argument("--rank", type=int, default=2)
a = p.parse_args()

x = np.linspace(-10, 10, 400)
t = np.linspace(0, 4 * np.pi, 200)
dt = t[1] - t[0]
Xg, T = np.meshgrid(x, t)
f1 = 1 / np.cosh(Xg + 3 - (0 if a.standing else T)) * np.exp(2.3j * T)
f2 = (np.tanh(Xg) / np.cosh(Xg)) * 2 * np.exp(2.8j * T)
X = (f1 + f2).T  # space x time

dmd = DMD(X, rank=a.rank, dt=dt, truncation=None)
X1, X2 = X[:, :-1], X[:, 1:]
U, S, Vh = np.linalg.svd(X1, full_matrices=False)
r = a.rank
Ur, Sr, Vr = U[:, :r], S[:r], Vh[:r].conj().T
lam, W = np.linalg.eig(Ur.conj().T @ X2 @ Vr / Sr)
Phi = X2 @ Vr / Sr @ W
omega = np.log(lam) / dt
b = np.linalg.lstsq(Phi, X[:, 0], rcond=None)[0]
X_hand = Phi @ (b[:, None] * np.exp(np.outer(omega, t)))
print("library omega:", np.sort_complex(dmd.omega))
print("hand omega:   ", np.sort_complex(omega))
print("relative reconstruction error: library", np.linalg.norm(dmd.reconstruct() - X) / np.linalg.norm(X),
      " hand", np.linalg.norm(X_hand - X) / np.linalg.norm(X))

fig = plt.figure(figsize=(12, 8))
for k, (data, title) in enumerate(((f1.T, "f1"), (f2.T, "f2"), (X, "f"), (dmd.reconstruct(), f"DMD rank {r}"))):
    ax = fig.add_subplot(2, 2, k + 1, projection="3d")
    ax.plot_surface(Xg, T, data.real.T, cmap="gray", linewidth=0, rstride=4, cstride=8)
    ax.set(title=title, xlabel="x", ylabel="t")
    ax.view_init(60, -20)
fig, axes = plt.subplots(1, 3, figsize=(13, 4))
axes[0].semilogy(dmd.singular_values / dmd.singular_values.sum(), "r+")
axes[0].set_title("normalised singular values")
axes[1].plot(dmd.omega.real, dmd.omega.imag, "+r", omega.real, omega.imag, "ob", mfc="none")
axes[1].set_title("continuous-time eigenvalues")
axes[2].plot(x, dmd.modes.real)
axes[2].set_title("DMD modes (real part)")
fig, ax = plt.subplots()
ax.plot(t, dmd.time_dynamics(t).real.T)
ax.set(title="time dynamics", xlabel="t")
plt.show()
