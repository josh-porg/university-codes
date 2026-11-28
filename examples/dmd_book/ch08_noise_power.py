"""DMD book chapter 8 (``DMD_eig``, ``FFTDMD_spectrum``, ``LDS_DMD_eig``, ``SVHT_cylinder``,
``optimal_SVHT_coef``): noise, power spectra and truncation.

1. FFT versus DMD spectrum of two noisy sines (7 and 13 Hz) after
   shift-stacking the signal 500 times; DMD power is ``2 |b| / sqrt(s)``.
2. A 2 x 2 linear system with noise (sigma 0.5): eigenvalues from an
   ensemble of 500 noisy realisations with exact, forward-backward and
   total-least-squares DMD (:func:`unicodes.decomposition.dmd_eigenvalues`);
   exact DMD is biased towards the origin.
3. With ``CYLINDER_ALL.mat``: optimal hard threshold of the singular values
   of noisy cylinder data (:func:`unicodes.decomposition.optimal_svht_coef`),
   and of pure noise of the same size::

    python ch08_noise_power.py [CYLINDER_ALL.mat]
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from unicodes.decomposition import dmd_eigenvalues, optimal_svht_coef, time_delay_stack

p = argparse.ArgumentParser()
p.add_argument("cylinder", nargs="?")
p.add_argument("--ensemble", type=int, default=500)
a = p.parse_args()
rng = np.random.default_rng(0)

# 1. FFT vs DMD
dt = 0.01
t = np.arange(0, 10 + dt / 2, dt)
xclean = 14 * np.sin(7 * 2 * np.pi * t) + 5 * np.sin(13 * 2 * np.pi * t)
x = xclean + 10 * rng.standard_normal(t.size)
fig, ax = plt.subplots(2, 1, figsize=(7, 5))
ax[0].plot(t, x, color=(0.7, 0.3, 0.3), lw=0.9, label="Noisy signal")
ax[0].plot(t, xclean, "k", lw=1.5, label="Clean signal")
ax[0].legend()
N = t.size
xhat = np.fft.fft(x)
freqs = np.arange(N // 2 + 1) / (N * dt)
ax[1].plot(freqs, np.abs(xhat[:N // 2 + 1]) * 2 / N, "k", lw=1.2, label="FFT")
s = 500
X = time_delay_stack(x, s)
U, S, Vh = np.linalg.svd(X[:, :-1], full_matrices=False)
r = 50
At = U[:, :r].T @ X[:, 1:] @ Vh[:r].T / S[:r]
lam, W = np.linalg.eig(At)
Phi = X[:, 1:] @ Vh[:r].T / S[:r] @ W
b = np.linalg.lstsq(Phi, X[:, 0], rcond=None)[0]
dmd_f = np.abs(np.log(lam).imag / dt / (2 * np.pi))
power = np.abs(b) * 2 / np.sqrt(s)
ax[1].scatter(dmd_f, power, edgecolors="r", facecolors="none", label="DMD")
ax[1].legend()
top = np.argsort(power)[::-1][:4]
print("strongest DMD frequencies (Hz):", np.round(dmd_f[top], 3), "power", np.round(power[top], 2))

# 2. noisy linear system
A = np.array([[1, 1], [-1, 2]]) / np.sqrt(3)
m = 100
Xl = np.zeros((2, m))
Xl[:, 0] = [0.5, 1]
for k in range(1, m):
    Xl[:, k] = A @ Xl[:, k - 1]
true = np.linalg.eigvals(A)
ens = {f: [] for f in ("dmd", "fbdmd", "tlsdmd")}
for _ in range(a.ensemble):
    Y = Xl + 0.5 * rng.standard_normal(Xl.shape)
    for f in ens:
        ens[f].append(dmd_eigenvalues(Y, f))
plt.figure()
th = np.linspace(0, 2 * np.pi, 200)
plt.plot(np.cos(th), np.sin(th), "k--")
for f, style in (("dmd", "b."), ("tlsdmd", "o"), ("fbdmd", "k+")):
    lam = np.concatenate(ens[f])
    upper = lam[lam.imag > 0]
    plt.plot(lam.real, lam.imag, style, ms=3, label=f, mfc="none")
    print(f"{f:7s}: mean of upper eigenvalue {upper.mean():.4f} (true {true[true.imag > 0][0]:.4f})")
plt.plot(true.real, true.imag, "r^", label="true")
plt.xlim(true.real[0] - 0.45, true.real[0] + 0.25)
plt.ylim(0, 0.8)
plt.gca().set_aspect("equal")
plt.legend()
plt.title("Noise = 0.5")

# 3. SVHT
if a.cylinder:
    from cylinder_plot import load_cylinder, plot_cylinder

    VORTALL, _, nx, ny = load_cylinder(a.cylinder)
    Y = VORTALL + 0.01 * rng.standard_normal(VORTALL.shape)
    U, sig, _ = np.linalg.svd(Y, full_matrices=False)
    tau = optimal_svht_coef(Y.shape[1] / Y.shape[0]) * np.median(sig)
    print(f"SVHT keeps {np.sum(sig > tau)} of {sig.size} singular values")
    plt.figure()
    plt.semilogy(sig, "k-o", lw=1.2)
    plt.semilogy(np.flatnonzero(sig > tau), sig[sig > tau], "ro")
    plt.axis([0, 150, 1, 1e4])
    plot_cylinder(U[:, 26] / np.abs(U[:, 26]).max() * 5, "singular vector 27")
    noise = 0.01 * rng.standard_normal(VORTALL.shape)
    sn = np.linalg.svd(noise, compute_uv=False)
    thresh = optimal_svht_coef(noise.shape[1] / noise.shape[0]) * np.median(sn)
    print(f"pure noise: SVHT keeps {np.sum(sn > thresh)} singular values")
    fig, ax = plt.subplots(1, 2)
    ax[0].semilogy(sn, "b")
    ax[1].plot(np.cumsum(sn / sn.sum()))
else:
    print("SVHT part skipped (pass CYLINDER_ALL.mat)")
plt.show()
