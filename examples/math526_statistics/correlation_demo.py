"""MATH 526: correlation coefficients of noisy positively, negatively and un-correlated signals
(``correlation_demonstration``, ``plotSingalVTime``, ``saveAllFig``).

``awgn(x, snr_dB)`` (Communications Toolbox) is reproduced: white Gaussian
noise scaled to the given signal-to-noise ratio of each signal. For the
all-zero signal ``awgn`` adds unit-power noise relative to 0 dBW; the same
convention is used. ``--save DIR`` writes every figure as PNG (``saveAllFig``).
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("--save", type=Path, default=None)
a = p.parse_args()
rng = np.random.default_rng(0)


def awgn(x, snr_db):
    x = np.atleast_2d(x)
    power = np.mean(x**2, axis=1, keepdims=True)
    power = np.where(power > 0, power, 1.0)
    return x + np.sqrt(power / 10 ** (snr_db / 10)) * rng.standard_normal(x.shape)


t = np.arange(1, 1001)
x = np.sin(t * 3 / 17)
ideal = np.vstack([x, -x, np.zeros_like(x)])


def correlations(snr_db):
    noisy = awgn(ideal, snr_db)
    r = lambda u, v: np.corrcoef(u, v)[0, 1]  # noqa: E731
    return noisy, [r(ideal[0], noisy[0]), r(ideal[0], noisy[1]), r(ideal[0], noisy[2]), r(noisy[0], noisy[1]),
                   r(noisy[0], noisy[2])]


print("correlation with t:", [np.corrcoef(s, t)[0, 1] if s.any() else np.nan for s in ideal])
noisy, rho = correlations(20)
print("rho at 20 dB (pos, neg, uncorrelated, noisy-neg, noisy-uncorrelated):", np.round(rho, 3))
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for ax, k, title in zip(axes, range(3), ("positively correlated", "negatively correlated", "uncorrelated")):
    ax.plot(ideal[0], noisy[k], "x", ms=3)
    ax.set(title=title, xlabel="sin(3t/17)", aspect="equal")
fig, ax = plt.subplots()
ax.plot(t, noisy.T, "-x", ms=2)
ax.set(title="Noisy data", xlabel="Time", ylabel="Amplitude")

levels = np.linspace(0, 50, 200)
rhos = np.array([correlations(s)[1] for s in levels])
fig, ax = plt.subplots()
ax.plot(levels, rhos, "-x", ms=3)
ax.legend(["pos", "neg", "uncorrelated", "noisy pos/neg", "noisy pos/uncorrelated"])
ax.set(title="Correlation convergence", xlabel="Signal-to-noise ratio (dB)", ylabel="Correlation coefficient")
if a.save:
    a.save.mkdir(parents=True, exist_ok=True)
    for n in plt.get_fignums():
        plt.figure(n).savefig(a.save / f"figure_{n}.png")
plt.show()
