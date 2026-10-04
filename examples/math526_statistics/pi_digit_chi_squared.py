"""MATH 526: chi-squared test of the digit distribution of pi and e
(``didgital_distribution_chi_squared``, ``didgital_distribution_chi_squared_graphs``).

Digits come from ``mpmath`` (installed with sympy). The MATLAB went to
10^8 digits; the default stops at 10^5 (``--max-digits`` to change). Note
the MATLAB ``chiAndGraph`` always used pi, also for the "e" case; both
constants are tested here.
"""

import argparse

import matplotlib.pyplot as plt
import mpmath
import numpy as np
from scipy.stats import chi2

p = argparse.ArgumentParser()
p.add_argument("--max-digits", type=int, default=100_000)
a = p.parse_args()


def digits(constant, n):
    mpmath.mp.dps = n + 10
    s = mpmath.nstr(constant(), n + 5, strip_zeros=False).replace(".", "")
    return np.frombuffer(s[:n].encode(), dtype=np.uint8) - ord("0")


def chi_squared(d):
    counts = np.bincount(d, minlength=10)
    expected = len(d) / 10
    return np.sum((counts - expected) ** 2 / expected), counts / len(d)


all_pi = digits(lambda: mpmath.pi, a.max_digits)
ns = np.unique(np.floor(np.logspace(1, np.log10(a.max_digits), 40)).astype(int))
stats = np.array([chi_squared(all_pi[:n])[0] for n in ns])
dists = np.array([chi_squared(all_pi[:n])[1] for n in ns])
for name, const in (("pi", lambda: mpmath.pi), ("e", lambda: mpmath.e)):
    c, dist = chi_squared(digits(const, 3000))
    print(f"{name}, 3000 digits: chi^2 = {c:.3f}, p-value = {chi2.sf(c, 9):.3f}, frequencies {np.round(dist, 3)}")

fig, ax = plt.subplots()
ax.semilogx(ns, stats, "xr")
ax.axhline(chi2.ppf(0.95, 9), color="k", ls="--", label="95 % critical value")
ax.set(title=r"Convergence of $\chi^2$ for the digits of pi", xlabel="n", ylabel=r"$\chi^2$")
ax.legend()
fig, ax = plt.subplots()
ax.semilogx(ns, dists, "x")
ax.set(title="Convergence of digit frequencies", xlabel="Number of digits", ylabel="Frequency")
fig, ax = plt.subplots()
ax.hist(stats, 10)
ax.set(title=r"$\chi^2$ values", xlabel=r"$\chi^2$")
plt.show()
