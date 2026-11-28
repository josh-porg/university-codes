"""DMD book chapter 10 (``Algorithm_10_1`` - ``Algorithm_10_6``, ``Algorithm_10_7through11``,
``dmd_soliton_rhs``): Koopman observables for the nonlinear Schrodinger equation.

``i u_t + u_xx / 2 + |u|^2 u = 0`` is solved spectrally on a periodic
domain of length 30 (``u_t^ = -i k^2 u^ / 2 + i FFT(|u|^2 u)``). Rank-10 DMD
of the state ``u`` alone is compared with DMD on the observables
``g1 = (u, |u|^2 u)`` (the nonlinearity of the equation, which gives a
nearly exact Koopman embedding) and ``g2 = (u, |u|^2)``. Then the N = 2
breather: DMD, extended DMD (``X2 pinv(X1)`` and ``X2 X1^T pinv(X1 X1^T)``)
and kernel DMD with four kernels (:func:`unicodes.decomposition.kernel_dmd`).
The breather is re-solved with 2000 (extended DMD) and 200 (kernel DMD)
snapshots as in the book (``--edmd-slices``, ``--kernel-slices``).
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from unicodes.decomposition import kernel_dmd

p = argparse.ArgumentParser()
p.add_argument("--kernel-slices", type=int, default=200)
p.add_argument("--edmd-slices", type=int, default=2000)
a = p.parse_args()


def nls(n, L, t, N):
    x = np.linspace(-L / 2, L / 2, n + 1)[:n]
    k = 2 * np.pi / L * np.r_[np.arange(n // 2), np.arange(-n // 2, 0)]

    def rhs(_, ut):
        u = np.fft.ifft(ut)
        return -0.5j * k**2 * ut + 1j * np.fft.fft(np.abs(u) ** 2 * u)

    u0 = N / np.cosh(x)
    sol = solve_ivp(rhs, (t[0], t[-1]), np.fft.fft(u0).astype(complex), t_eval=t, method="DOP853", rtol=1e-10, atol=1e-10)
    return x, np.fft.ifft(sol.y, axis=0), u0  # (n, len(t))


def dmd_fit(Y1, Y2, r, y0, t, dt):
    U, S, Vh = np.linalg.svd(Y1, full_matrices=False)
    At = U[:, :r].conj().T @ Y2 @ Vh[:r].conj().T / S[:r]
    lam, W = np.linalg.eig(At)
    Phi = Y2 @ Vh[:r].conj().T / S[:r] @ W
    omega = np.log(lam) / dt
    b = np.linalg.lstsq(Phi, y0, rcond=None)[0]
    return Phi, omega, Phi @ (b[:, None] * np.exp(np.outer(omega, t))), S


# Algorithms 10.1-10.6: N = 2 sech initial condition, t in [0, pi], 21 slices
n, L = 512, 30
t = np.linspace(0, np.pi, 21)
dt = t[1] - t[0]
x, X, q = nls(n, L, t, 2)
X1, X2 = X[:, :-1], X[:, 1:]
errs = {}
_, _, q_dmd, _ = dmd_fit(X1, X2, 10, X[:, 0], t, dt)
errs["DMD on u"] = q_dmd
_, _, q1, _ = dmd_fit(np.vstack([X1, X1 * np.abs(X1) ** 2]), np.vstack([X2, X2 * np.abs(X2) ** 2]), 10,
                      np.r_[q, q * np.abs(q) ** 2], t, dt)
errs["g1 = (u, |u|^2 u)"] = q1[:n]
_, _, q2, _ = dmd_fit(np.vstack([X1, np.abs(X1) ** 2]), np.vstack([X2, np.abs(X2) ** 2]), 10, np.r_[q, np.abs(q) ** 2], t, dt)
errs["g2 = (u, |u|^2)"] = q2[:n]
fig = plt.figure(figsize=(12, 3.5))
for k, (name, data) in enumerate([("PDE", X)] + list(errs.items())):
    if name != "PDE":
        print(f"{name:18s}: relative error {np.linalg.norm(data - X) / np.linalg.norm(X):.3e}")
    ax = fig.add_subplot(1, 4, k + 1, projection="3d")
    T_, X_ = np.meshgrid(t, x)
    ax.plot_surface(X_, T_, np.abs(data), cmap="gray", linewidth=0)
    ax.set(title=name, xlim=(-15, 15), zlim=(0, 4))

# Algorithms 10.7-10.11: breather over t in [0, 2 pi]
x, U40, u0 = nls(512, 30, np.linspace(0, 2 * np.pi, 41), 2)
t40 = np.linspace(0, 2 * np.pi, 41)
Phi, omega, u_dmd, S = dmd_fit(U40[:, :-1], U40[:, 1:], 10, u0, t40, t40[1] - t40[0])
err = np.linalg.norm(u_dmd - U40, axis=0)
print(f"breather DMD (rank 10): error grows from {err[0]:.2e} to {err[-1]:.2e}")
fig, ax = plt.subplots(1, 3, figsize=(12, 3.5))
ax[0].plot(S, "ko")
ax[0].set_title("singular values")
ax[1].plot(t40, err)
ax[1].set_title("DMD error vs time")
ax[2].plot(x, np.abs(Phi[:, :4]))
ax[2].set(title="|DMD modes|", xlim=(-20, 20))

te = np.linspace(0, 2 * np.pi, a.edmd_slices + 1)
dte = te[1] - te[0]
_, Ue, _ = nls(256, 30, te, 2)
X1, X2 = Ue[:, :-1], Ue[:, 1:]
lam1 = np.linalg.eigvals(X2 @ np.linalg.pinv(X1))
lam2 = np.linalg.eigvals(X2 @ X1.T @ np.linalg.pinv(X1 @ X1.T))
fig, ax = plt.subplots(2, 1, num="extended and kernel DMD")
for axk in ax:
    axk.plot((np.log(lam1) / dte).real, (np.log(lam1) / dte).imag, "ro", label="X2 pinv(X1)")
    axk.plot((np.log(lam2) / dte).real, (np.log(lam2) / dte).imag, "ko", label="X2 X1^T pinv(X1 X1^T)")
ax[1].axis([-1000, 500, -1000, 1000])

tk = np.linspace(0, 2 * np.pi, a.kernel_slices + 1)
dtk = tk[1] - tk[0]
_, Uk, _ = nls(256, 30, tk, 2)
X1, X2 = Uk[:, :-1], Uk[:, 1:]
pw = 20
kernels = {
    "(1 + x^T y)^20": lambda u, v: (1 + u @ v) ** pw,
    "(1 + |x|^T |y|)^20": lambda u, v: (1 + np.abs(u) @ np.abs(v)) ** pw,
    "exp(-|(x-y)^T (x-y)|^2)": lambda u, v: np.exp(-np.abs((u - v) @ (u - v)) ** 2),
    "exp(-|x^T y|)": lambda u, v: np.exp(-np.abs(u @ v)),
}
fig, axs = plt.subplots(2, 2, figsize=(10, 7))
for axk, (name, kern) in zip(axs.flat, kernels.items()):
    with np.errstate(over="ignore", invalid="ignore"):
        lam = kernel_dmd(X1, X2, kern, tol=1e-10)
    om = np.log(lam.astype(complex)) / dtk
    axk.plot(om.real, om.imag, "ko")
    axk.set(title=name, xlim=(-3000, 3000), ylim=(-100, 100))
    axk.grid(True)
    print(f"kernel {name:24s}: {lam.size} eigenvalues")
gauss = kernel_dmd(X1, X2, kernels["exp(-|(x-y)^T (x-y)|^2)"], tol=1e-10)
om = np.log(gauss.astype(complex)) / dtk
for axk in ax:
    axk.plot(om.real, om.imag, "bo", label="kernel DMD")
ax[0].legend()
plt.show()
