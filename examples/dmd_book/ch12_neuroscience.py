"""DMD book chapter 12 (``DMD_ECoG``): DMD of electrocorticography recordings.

Needs ``ecog_window.mat`` (the book's ``DATA/NEURO``; ``X`` with 59
electrodes, ``dt``). The data are shift-stacked 17 times and rank-200 DMD
gives the eigenvalues and a mode spectrum (with the alternative mode
scaling ``A_hat = S^(-1/2) A_tilde S^(1/2)``), compared with the FFT power
spectra of the electrodes. The book averaged ``fftp(c, :)`` over the last
electrode only; the mean over all electrodes is drawn here::

    python ch12_neuroscience.py ecog_window.mat
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from unicodes.decomposition import time_delay_stack
from unicodes.io import load_mat

p = argparse.ArgumentParser()
p.add_argument("data")
p.add_argument("--rank", type=int, default=200)
p.add_argument("--stacks", type=int, default=17)
a = p.parse_args()

d = load_mat(a.data)
X0 = np.asarray(d["X"], float)
dt = float(np.squeeze(d["dt"]))
Xaug = time_delay_stack(X0, a.stacks)
X, Y = Xaug[:, :-1], Xaug[:, 1:]
U, S, Vh = np.linalg.svd(X, full_matrices=False)
r = a.rank
Ur, Sr, Vr = U[:, :r], S[:r], Vh[:r].T
At = Ur.T @ Y @ Vr / Sr
lam = np.linalg.eigvals(At)
omega = np.log(lam.astype(complex)) / dt / (2 * np.pi)

fig, ax = plt.subplots(1, 2, figsize=(10, 4.5))
ax[0].plot(lam.real, lam.imag, "k.")
th = np.linspace(0, 2 * np.pi, 200)
ax[0].plot(np.cos(th), np.sin(th), "k--")
ax[0].set(aspect="equal", xlim=(-1.2, 1.2), ylim=(-1.2, 1.2))
ax[1].plot(omega.real, omega.imag, "k.")
ax[1].axvline(0, color="k", ls="--")
ax[1].set(xlim=(-8, 2), ylim=(-170, 170))

s_half = np.sqrt(Sr)
lam_h, W_h = np.linalg.eig((At * s_half[None, :]) / s_half[:, None])  # S^-1/2 A_tilde S^1/2
W = s_half[:, None] * W_h
Phi = Y @ Vr / Sr @ W
omega_h = np.log(lam_h.astype(complex)) / dt / (2 * np.pi)
P = np.real(np.sum(np.conj(Phi) * Phi, axis=0))
fig, ax = plt.subplots(1, 2, figsize=(10, 4.5))
ax[0].stem(np.abs(omega_h.imag), P)
ax[0].set(xlim=(0, 150), title="DMD spectrum")
nfft = 2 ** int(np.ceil(np.log2(X.shape[1])))
f = 1 / dt / 2 * np.linspace(0, 1, nfft // 2 + 1)
fftp = np.fft.fft(X0[:, :X.shape[1]], nfft, axis=1)
for c in range(X0.shape[0]):
    ax[1].plot(f, 2 * np.abs(fftp[c, :nfft // 2 + 1]), color=(0.6, 0.6, 0.6))
ax[1].plot(f, 2 * np.abs(fftp[:, :nfft // 2 + 1]).mean(axis=0), "k", lw=2)
ax[1].set(xlim=(0, 150), ylim=(0, 400), title="FFT power")
plt.show()
