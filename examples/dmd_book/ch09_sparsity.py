"""DMD book chapter 9 (``Algorithm_9_1``, ``cosamp``, ``EX1_TORUS``: ``runExample``, ``getParms``,
``getSparseData``, ``plotData``, ``computeFFTModes``, ``computePODModes``, ``computeDMDModes``,
``projectData``, ``compressedDMD``; ``EX2_CYLINDER``: ``cDMD_p40_r21_VORT``, ``csDMD_p1000_r21_VORT``,
``plotCylinderNoSave``): compressed sensing and compressed DMD.

1. Two sines (73 and 531 Hz) sampled at 256 random points of 4096 are
   recovered by CoSaMP in the DCT basis (:func:`unicodes.decomposition.cosamp`).
   The sines are only nearly sparse in the DCT (they leak into neighbouring
   coefficients), so a 10-sparse fit leaves a 20-30 % error. The book fitted
   ``Psi(perm, :) s = y`` with the forward DCT matrix and reconstructed with the
   inverse transform; one basis (``x = Psi s``) is used for both here.
2. Torus example: K = 5 damped oscillating Fourier modes on a 128 x 128
   grid mapped onto a torus. FFT, POD and DMD modes of the full data, then
   compressed DMD from M = 15 Gaussian random projections: the eigenvalues
   come from the projected data and the full modes from the full data with
   the projected-data DMD vectors (``PhiXr``), or from the compressed modes
   alone by CoSaMP in the Fourier basis. The book drew the Fourier indices
   from 1..9 (MATLAB), which includes the mean; a real mean mode oscillating
   in time is a standing wave of rank one that DMD cannot represent (see
   ``ch07_time_delay.py``), so the indices start at the first wavenumber here.
3. With ``--cylinder CYLINDER_ALL.mat``: compressed DMD of the cylinder wake
   from ``--p`` random projections (40 in ``cDMD``; ``csDMD`` used 1000 and
   reconstructed the modes by compressed sensing, ``--cs``).

The book's ``freezeColors`` (several colormaps per MATLAB figure) is not
needed with matplotlib::

    python ch09_sparsity.py
    python ch09_sparsity.py --cylinder CYLINDER_ALL.mat --p 1000 --cs
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import dct

from unicodes.decomposition import cosamp

p = argparse.ArgumentParser()
p.add_argument("--cylinder")
p.add_argument("--p", type=int, default=40)
p.add_argument("--cs", action="store_true", help="reconstruct cylinder modes by compressed sensing")
p.add_argument("--seed", type=int, default=1)
a = p.parse_args()
rng = np.random.default_rng(a.seed)

# 1. Algorithm 9.1
n, p_ = 4096, 256
t = np.linspace(0, 1, n)
x = np.sin(73 * 2 * np.pi * t) + np.sin(531 * 2 * np.pi * t)
perm = rng.choice(n, p_, replace=False)
Psi = dct(np.eye(n), norm="ortho", axis=0).T  # columns are DCT basis vectors: x = Psi s
s = cosamp(Psi[perm], x[perm], 10, 1e-10, 10)
x_rec = Psi @ s
print(f"9.1: CoSaMP reconstruction error {np.linalg.norm(x_rec - x) / np.linalg.norm(x):.3e} from {p_} of {n} samples")
plt.figure()
plt.plot(t, x, t, x_rec, "--")
plt.xlim(0.4, 0.5)

# 2. torus
n, K, T, dt, r, M = 128, 5, 2.0, 0.01, 10, 15
xt = np.zeros((n, n))
I, J, IC, F, damping = [], [], [], [], []
for _ in range(K):
    while True:
        i, j = rng.integers(1, n // 15 + 1, 2)  # not the mean (0, 0): see the docstring
        if xt[i, j] == 0:
            break
    I.append(i)
    J.append(j)
    IC.append(rng.standard_normal())
    xt[i, j] = IC[-1]
    F.append(np.sqrt(4 * rng.random()))
    damping.append(-rng.random() * 0.1)
I, J, IC, F, damping = map(np.array, (I, J, IC, F, damping))
times = np.arange(0, T + dt / 2, dt)
XT = np.zeros((n * n, times.size), complex)
for c, tt in enumerate(times):
    xt = np.zeros((n, n), complex)
    xt[I, J] = np.exp(damping * tt) * np.exp(1j * 2 * np.pi * F * tt) * IC
    XT[:, c] = xt.ravel()
XDAT = np.real(np.fft.ifft2(XT.T.reshape(-1, n, n))).reshape(times.size, -1).T

T1, T2 = np.meshgrid(np.linspace(0, 2 * np.pi, n), np.linspace(0, 2 * np.pi, n))
R = 2 + np.cos(T2)
TX, TY, TZ = np.cos(T1) * R, np.sin(T1) * R, np.sin(T2)


def torus(ax, field, title=""):
    f = np.real(field.reshape(n, n))
    ax.plot_surface(TX, TY, TZ, facecolors=plt.cm.jet((f - f.min()) / (np.ptp(f) + 1e-30)), linewidth=0,
                    rstride=2, cstride=2)
    ax.view_init(32, 10)
    ax.set_axis_off()
    ax.set_title(title, fontsize=8)


fig = plt.figure(figsize=(12, 4))
for k in range(K):  # FFT modes
    for row in range(2):
        x2t = np.zeros((n, n), complex)
        x2t[I[k], J[k]] = 1j**row
        torus(fig.add_subplot(2, K, row * K + k + 1, projection="3d"), np.fft.ifft2(x2t), f"FFT mode {k + 1}")

U, S, Vh = np.linalg.svd(XDAT, full_matrices=False)
print(f"POD: {np.sum(np.cumsum(S) / S.sum() < 0.999) + 1} modes hold 99.9 % of the singular-value sum (2K = {2 * K})")


def dmd(X, Xref=None):
    X1, X2 = X[:, :-1], X[:, 1:]
    U, S, Vh = np.linalg.svd(X1, full_matrices=False)
    At = U[:, :r].conj().T @ X2 @ Vh[:r].conj().T / S[:r]
    lam, W = np.linalg.eig(At)
    Phi = X2 @ Vh[:r].conj().T / S[:r] @ W
    PhiX = None if Xref is None else Xref[:, 1:] @ Vh[:r].conj().T / S[:r] @ W
    return np.log(lam) / dt, Phi, PhiX


true = np.concatenate([damping + 2j * np.pi * F, damping - 2j * np.pi * F])
eig_full, Phi_full, _ = dmd(XDAT)
C = rng.standard_normal((M, n * n))
Theta = np.array([np.fft.ifft2(row.reshape(n, n)).ravel() for row in C])
YDAT = C @ XDAT
eig_c, PhiY, PhiXr = dmd(YDAT, XDAT)


def match_err(lam):
    return np.max(np.min(np.abs(lam[:, None] - true[None, :]), axis=0))


print(f"torus: max eigenvalue error, full DMD {match_err(eig_full):.2e}, compressed DMD ({M} projections) {match_err(eig_c):.2e}")
plt.figure()
plt.scatter(true.real, true.imag, marker="x", c="k", label="True eigenvalues")
plt.scatter(eig_full.real, eig_full.imag, facecolors="none", edgecolors="b", label="Exact DMD eigenvalues")
plt.scatter(eig_c.real, eig_c.imag, marker="d", facecolors="none", edgecolors="r", label="Compressed DMD eigenvalues")
plt.legend()
# compressed sensing reconstruction of the modes from PhiY
Phi_cs = []
for i in range(0, 2 * K, 2):
    ahat = cosamp(Theta, PhiY[:, i], 2, 1e-5, 100)
    m_ = np.fft.ifft2(ahat.reshape(n, n))
    Phi_cs += [m_.real.ravel(), m_.imag.ravel()]
fig = plt.figure(figsize=(12, 4))
for k in range(K):
    torus(fig.add_subplot(2, K, k + 1, projection="3d"), PhiXr[:, 2 * k] + PhiXr[:, 2 * k + 1], "compressed DMD")
    torus(fig.add_subplot(2, K, K + k + 1, projection="3d"), Phi_cs[2 * k], "CS-DMD")

# 3. cylinder
if a.cylinder:
    from cylinder_plot import load_cylinder, plot_cylinder

    VORTALL, VORTEXTRA, nx, ny = load_cylinder(a.cylinder)
    X = VORTALL
    X2 = np.hstack([VORTALL[:, 1:], VORTEXTRA]) if VORTEXTRA.size else VORTALL[:, 1:]
    X = X[:, :X2.shape[1]]
    C = rng.standard_normal((a.p, X.shape[0]))
    Y, Y2 = C @ X, C @ X2
    S = np.linalg.svd(X, compute_uv=False)
    UY, SY, VYh = np.linalg.svd(Y, full_matrices=False)
    plt.figure()
    plt.semilogy(S[:20], "-ok", SY[:20], "or")
    plt.title("singular values: full (black), projected (red)")
    r = 21
    At = UY[:, :r].T @ Y2 @ VYh[:r].T / SY[:r]
    lam, WY = np.linalg.eig(At)
    PhiY = Y2 @ VYh[:r].T / SY[:r] @ WY
    PhiXt = X2 @ VYh[:r].T / SY[:r] @ WY
    for i in range(9, r, 2):
        plot_cylinder(PhiXt[:, i].real / np.abs(PhiXt[:, i]).max() * 5, f"compressed DMD mode {i + 1}")
    if a.cs:
        Theta = np.array([np.fft.ifft2(row.reshape(ny, nx).T).T.ravel() for row in C])
        for i in range(9, r):
            for k in (25, 50, 75):
                ahat = cosamp(Theta, PhiY[:, i], k, 1e-5, 10)
                m_ = np.fft.ifft2(ahat.reshape(ny, nx).T)
                plot_cylinder(m_.real / np.abs(m_).max() * 5, f"CS mode {i + 1}, K = {k}")
plt.show()
