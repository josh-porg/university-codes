"""DMD book chapter 7 (``DMD_standingwave``, ``ERA``, ``ERA_test01``, ``HMM_DMD``): time-delay
coordinates, the eigensystem realization algorithm and hidden Markov models.

1. A single measurement of a standing wave ``sin t`` has rank one: DMD
   returns one real eigenvalue and cannot oscillate; one delay coordinate
   recovers ``omega = +-i``.
2. ERA (:func:`unicodes.decomposition.era`) on the impulse response
   ``e^{-0.01 t} + e^{-10 t}`` recovers both decay rates.
3. A weather Markov chain (sun/rain/clouds) with a dog's mood as the
   emission: ten-day forecast, simulated sequences, Viterbi decoding,
   counting estimates, posterior decoding and Viterbi training (the MATLAB
   Statistics Toolbox ``hmm*`` functions are written out here), the
   observability rank, and DMD of one-hot state data, whose eigenvalues
   are those of the transition matrix.
"""

import matplotlib.pyplot as plt
import numpy as np

from unicodes.decomposition import era

rng = np.random.default_rng(0)

# 1. standing wave
dt = 0.01
t = np.arange(0, 10 + dt / 2, dt)
x = np.sin(t)
X, X2 = x[None, :-1], x[None, 1:]
U, S, Vh = np.linalg.svd(X, full_matrices=False)
omega1 = np.log((U.T @ X2 @ Vh.T / S).item()) / dt
print(f"DMD on x alone: omega = {omega1:.4f} (real: no oscillation)")
basis = np.exp(omega1 * t)
xdmd = basis * (basis @ x) / (basis @ basis)  # least-squares amplitude over all samples
Xa, Xa2 = np.vstack([x[:-2], x[1:-1]]), np.vstack([x[1:-1], x[2:]])
U, S, Vh = np.linalg.svd(Xa, full_matrices=False)
lam, W = np.linalg.eig(U.T @ Xa2 @ Vh.T / S)
Om = np.log(lam) / dt
Phi = Xa2 @ Vh.T / S @ W
b = np.linalg.solve(Phi, Xa[:, 0])
xaug = (Phi @ (b[:, None] * np.exp(np.outer(Om, t))))[0].real
print(f"with one delay: omega = {np.round(Om, 6)}")
plt.figure(figsize=(6, 2.5))
plt.plot(t, x, label="Data")
plt.plot(t, xdmd, "r--", label="DMD")
plt.plot(t, xaug, "b--", lw=1.5, label="Augmented DMD")
plt.ylim(-1, 1)
plt.legend()

# 2. ERA
H = np.array([np.exp(-t[k:k + 5]) + np.exp(-3 * t[k:k + 5]) for k in range(5)])
print(f"rank of the 5 x 5 Hankel matrix of e^-t + e^-3t: {np.linalg.matrix_rank(H)}")
t = np.arange(0, 100 + dt / 2, dt)
x = np.exp(-0.01 * t) + np.exp(-10 * t)
r = np.linalg.matrix_rank(np.vstack([x[:-2], x[1:-1]]))
A, B, C, D, hsv = era(x[None, None, :], 100, 100, r)
y = np.zeros_like(t)
state = np.zeros(r)
for k in range(t.size):  # discrete impulse response
    u = 1.0 if k == 0 else 0.0
    y[k] = (C @ state).item() + D.item() * u
    state = A @ state + B.ravel() * u
print("ERA continuous eigenvalues:", np.round(np.log(np.linalg.eigvals(A)) / dt, 6))
fig, ax = plt.subplots(1, 2, figsize=(9, 3))
for axk, lim in zip(ax, (1, 100)):
    axk.plot(t, x, "k", t, y, "r--")
    axk.set(xlim=(0, lim), xlabel="Time")
ax[0].legend(["Data", "ERA Model"])


# 3. hidden Markov model
T = np.array([[0.00, 0.40, 0.60], [0.25, 0.50, 0.25], [0.25, 0.25, 0.50]])
E = np.array([[1.0, 0.0], [0.0, 1.0], [0.2, 0.8]])


def hmm_generate(N, T, E):
    states, seq = np.zeros(N, int), np.zeros(N, int)
    s = 0  # MATLAB starts in state 1 before the first step
    for i in range(N):
        s = rng.choice(len(T), p=T[s])
        states[i], seq[i] = s, rng.choice(E.shape[1], p=E[s])
    return seq, states


def hmm_viterbi(seq, T, E):
    with np.errstate(divide="ignore"):
        lT, lE = np.log(T), np.log(E)
    v = lT[0] + lE[:, seq[0]]
    back = np.zeros((len(seq), len(T)), int)
    for i in range(1, len(seq)):
        cand = v[:, None] + lT
        back[i] = cand.argmax(0)
        v = cand.max(0) + lE[:, seq[i]]
    path = [int(v.argmax())]
    for i in range(len(seq) - 1, 0, -1):
        path.append(back[i, path[-1]])
    return np.array(path[::-1])


def hmm_estimate(seq, states, n_states, n_symbols):
    Tc = np.zeros((n_states, n_states))
    np.add.at(Tc, (states[:-1], states[1:]), 1)
    Ec = np.zeros((n_states, n_symbols))
    np.add.at(Ec, (states, seq), 1)
    with np.errstate(invalid="ignore"):
        return np.nan_to_num(Tc / Tc.sum(1, keepdims=True)), np.nan_to_num(Ec / Ec.sum(1, keepdims=True))


def hmm_decode(seq, T, E):
    """Posterior state probabilities (scaled forward-backward)."""
    N, n = len(seq), len(T)
    f, c = np.zeros((N, n)), np.zeros(N)
    f[0] = T[0] * E[:, seq[0]]
    c[0] = f[0].sum()
    f[0] /= c[0]
    for i in range(1, N):
        f[i] = (f[i - 1] @ T) * E[:, seq[i]]
        c[i] = f[i].sum()
        f[i] /= c[i]
    bwd = np.ones((N, n))
    for i in range(N - 2, -1, -1):
        bwd[i] = T @ (E[:, seq[i + 1]] * bwd[i + 1]) / c[i + 1]
    return (f * bwd).T


def hmm_train_viterbi(seq, T, E, iterations=50):
    for _ in range(iterations):
        path = hmm_viterbi(seq, T, E)
        T_new, E_new = hmm_estimate(seq, path, len(T), E.shape[1])
        T_new = np.where(T_new.sum(1, keepdims=True) > 0, T_new, T)
        E_new = np.where(E_new.sum(1, keepdims=True) > 0, E_new, E)
        if np.allclose(T_new, T) and np.allclose(E_new, E):
            break
        T, E = T_new, E_new
    return T, E


w, V = np.linalg.eig(T.T)
pi_ss = np.real(V[:, np.argmax(w.real)])
print("stationary distribution [S R C]:", np.round(pi_ss / pi_ss.sum(), 4))
x_hist = [np.array([0.0, 1, 0])]
for _ in range(10):
    x_hist.append(x_hist[-1] @ T)
x_hist = np.array(x_hist)
fig, ax = plt.subplots(2, 1)
ax[0].plot(x_hist)
ax[0].legend(["Sunny", "Rainy", "Cloudy"])
ax[1].plot(x_hist @ E)
ax[1].legend(["Happy Dog", "Grumpy Dog"])
ax[1].set_xlabel("# of days from today")

seq, states = hmm_generate(1000, T, E)
likely = hmm_viterbi(seq, T, E)
print(f"Viterbi decoding correct for {np.mean(likely == states):.1%} of days")
T_est, E_est = hmm_estimate(seq, states, 3, 2)
print("estimated T:\n", np.round(T_est, 3), "\nestimated E:\n", np.round(E_est, 3))
pstates = hmm_decode(seq, T, E)
print(f"posterior decoding correct for {np.mean(pstates.argmax(0) == states):.1%} of days")
T2, E2 = hmm_train_viterbi(seq, T + 0.001 * rng.standard_normal(T.shape).clip(0), np.abs(E + 0.001 * rng.standard_normal(E.shape)))
print("Viterbi-trained T:\n", np.round(T2, 3))
plt.figure()
plt.plot(states, "k-")
plt.plot(likely, "rx")


def obsv_rank(A, C):
    return np.linalg.matrix_rank(np.vstack([C @ np.linalg.matrix_power(A, k) for k in range(len(A))]))


E_bad = np.array([[0.75, 0.25], [0.5, 0.5], [0.5, 0.5]])
print(f"observability rank: dog emissions {obsv_rank(T.T, E.T)}, ambiguous emissions {obsv_rank(T.T, E_bad.T)}")

# DMD of one-hot states (E = identity)
seq, _ = hmm_generate(10000, T, np.eye(3))
Xs = np.zeros((3, seq.size))
Xs[seq, np.arange(seq.size)] = 1
Y, Xs = Xs[:, 1:], Xs[:, :-1]
U, S, Vh = np.linalg.svd(Xs, full_matrices=False)
At = U.T @ Y @ Vh.T / S
print("DMD eigenvalues of one-hot data:", np.round(np.sort(np.linalg.eigvals(At).real), 3),
      " transition matrix:", np.round(np.sort(np.linalg.eigvals(T).real), 3))
plt.show()
