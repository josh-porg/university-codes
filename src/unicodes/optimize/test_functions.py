"""Optimisation benchmark functions (Math 796 HW 1 ``optimization_test_functions``).

The formulas are those of the Surjanovic & Bingham Virtual Library of
Simulation Experiments (sfu.ca/~ssurjano), which the MATLAB files were taken
from. Every function takes ``x`` with the design variables on the first axis,
so ``f(np.array([X, Y]))`` evaluates a mesh grid. ``CATEGORIES`` groups them
as the library does.
"""

from __future__ import annotations

import numpy as np


def _idx(x):
    return np.arange(1, x.shape[0] + 1).reshape((-1,) + (1,) * (x.ndim - 1))


def ackley(x, a=20.0, b=0.2, c=2 * np.pi):
    x = np.asarray(x, float)
    d = x.shape[0]
    return -a * np.exp(-b * np.sqrt((x**2).sum(0) / d)) - np.exp(np.cos(c * x).sum(0) / d) + a + np.e


def drop(x):
    r2 = x[0] ** 2 + x[1] ** 2
    return -(1 + np.cos(12 * np.sqrt(r2))) / (0.5 * r2 + 2)


def griewank(x):
    x = np.asarray(x, float)
    return (x**2).sum(0) / 4000 - np.prod(np.cos(x / np.sqrt(_idx(x))), axis=0) + 1


def langermann(x, c=(1, 2, 5, 2, 3), A=((3, 5), (5, 2), (2, 1), (1, 4), (7, 9))):
    x = np.asarray(x, float)
    out = 0.0
    for ci, Ai in zip(c, A):
        inner = sum((x[j] - Ai[j]) ** 2 for j in range(x.shape[0]))
        out = out + ci * np.exp(-inner / np.pi) * np.cos(np.pi * inner)
    return out


def levy(x):
    x = np.asarray(x, float)
    w = 1 + (x - 1) / 4
    s = ((w[:-1] - 1) ** 2 * (1 + 10 * np.sin(np.pi * w[:-1] + 1) ** 2)).sum(0)
    return np.sin(np.pi * w[0]) ** 2 + s + (w[-1] - 1) ** 2 * (1 + np.sin(2 * np.pi * w[-1]) ** 2)


def levy13(x):
    x1, x2 = x[0], x[1]
    return (np.sin(3 * np.pi * x1) ** 2 + (x1 - 1) ** 2 * (1 + np.sin(3 * np.pi * x2) ** 2)
            + (x2 - 1) ** 2 * (1 + np.sin(2 * np.pi * x2) ** 2))


def rastrigin(x):
    x = np.asarray(x, float)
    return 10 * x.shape[0] + (x**2 - 10 * np.cos(2 * np.pi * x)).sum(0)


def schaffer2(x):
    r2 = x[0] ** 2 + x[1] ** 2
    return 0.5 + (np.sin(x[0] ** 2 - x[1] ** 2) ** 2 - 0.5) / (1 + 0.001 * r2) ** 2


def shubert(x):
    i = np.arange(1, 6).reshape((-1,) + (1,) * (np.ndim(x[0])))
    return (i * np.cos((i + 1) * x[0] + i)).sum(0) * (i * np.cos((i + 1) * x[1] + i)).sum(0)


def bohachevsky1(x):
    return x[0] ** 2 + 2 * x[1] ** 2 - 0.3 * np.cos(3 * np.pi * x[0]) - 0.4 * np.cos(4 * np.pi * x[1]) + 0.7


def perm0db(x, b=10.0):
    x = np.asarray(x, float)
    j = _idx(x)
    return sum((((j + b) * (x**i - (1 / j) ** i)).sum(0)) ** 2 for i in range(1, x.shape[0] + 1))


def rotated_hyper_ellipsoid(x):
    x = np.asarray(x, float)
    return np.cumsum(x**2, axis=0).sum(0)


def sum_squares(x):
    x = np.asarray(x, float)
    return (_idx(x) * x**2).sum(0)


def trid(x):
    x = np.asarray(x, float)
    return ((x - 1) ** 2).sum(0) - (x[1:] * x[:-1]).sum(0)


def dejong5(x):
    a = np.array([-32, -16, 0, 16, 32])
    a1, a2 = np.tile(a, 5), np.repeat(a, 5)
    s = sum(1 / (i + 1 + (x[0] - a1[i]) ** 6 + (x[1] - a2[i]) ** 6) for i in range(25))
    return 1 / (0.002 + s)


def easom(x):
    return -np.cos(x[0]) * np.cos(x[1]) * np.exp(-((x[0] - np.pi) ** 2) - (x[1] - np.pi) ** 2)


def michalewicz(x, m=10):
    x = np.asarray(x, float)
    return -(np.sin(x) * np.sin(_idx(x) * x**2 / np.pi) ** (2 * m)).sum(0)


def booth(x):
    return (x[0] + 2 * x[1] - 7) ** 2 + (2 * x[0] + x[1] - 5) ** 2


def matyas(x):
    return 0.26 * (x[0] ** 2 + x[1] ** 2) - 0.48 * x[0] * x[1]


def mccormick(x):
    return np.sin(x[0] + x[1]) + (x[0] - x[1]) ** 2 - 1.5 * x[0] + 2.5 * x[1] + 1


def power_sum(x, b=None):
    x = np.asarray(x, float)
    d = x.shape[0]
    b = {2: (8, 18), 4: (8, 18, 44, 114)}[d] if b is None else b
    return sum(((x**i).sum(0) - b[i - 1]) ** 2 for i in range(1, d + 1))


def zakharov(x):
    x = np.asarray(x, float)
    s2 = (0.5 * _idx(x) * x).sum(0)
    return (x**2).sum(0) + s2**2 + s2**4


def camel3(x):
    return 2 * x[0] ** 2 - 1.05 * x[0] ** 4 + x[0] ** 6 / 6 + x[0] * x[1] + x[1] ** 2


def camel6(x):
    return (4 - 2.1 * x[0] ** 2 + x[0] ** 4 / 3) * x[0] ** 2 + x[0] * x[1] + (-4 + 4 * x[1] ** 2) * x[1] ** 2


def dixon_price(x):
    x = np.asarray(x, float)
    return (x[0] - 1) ** 2 + (_idx(x)[1:] * (2 * x[1:] ** 2 - x[:-1]) ** 2).sum(0)


def rosenbrock(x):
    x = np.asarray(x, float)
    return (100 * (x[1:] - x[:-1] ** 2) ** 2 + (x[:-1] - 1) ** 2).sum(0)


def beale(x):
    x1, x2 = x[0], x[1]
    return (1.5 - x1 + x1 * x2) ** 2 + (2.25 - x1 + x1 * x2**2) ** 2 + (2.625 - x1 + x1 * x2**3) ** 2


def branin(x, a=1.0, b=5.1 / (4 * np.pi**2), c=5 / np.pi, r=6.0, s=10.0, t=1 / (8 * np.pi)):
    return a * (x[1] - b * x[0] ** 2 + c * x[0] - r) ** 2 + s * (1 - t) * np.cos(x[0]) + s


def branin_modified(x):
    return branin(x) + 5 * x[0]


def goldstein_price(x):
    x1, x2 = x[0], x[1]
    f1 = 1 + (x1 + x2 + 1) ** 2 * (19 - 14 * x1 + 3 * x1**2 - 14 * x2 + 6 * x1 * x2 + 3 * x2**2)
    f2 = 30 + (2 * x1 - 3 * x2) ** 2 * (18 - 32 * x1 + 12 * x1**2 + 48 * x2 - 36 * x1 * x2 + 27 * x2**2)
    return f1 * f2


def permdb(x, b=0.5):
    x = np.asarray(x, float)
    j = _idx(x)
    return sum(((j**i + b) * ((x / j) ** i - 1)).sum(0) ** 2 for i in range(1, x.shape[0] + 1))


def styblinski_tang(x):
    x = np.asarray(x, float)
    return (x**4 - 16 * x**2 + 5 * x).sum(0) / 2


CATEGORIES = {
    "many_local_minima": [ackley, drop, griewank, langermann, levy, levy13, rastrigin, schaffer2, shubert],
    "bowl_shaped": [bohachevsky1, perm0db, rotated_hyper_ellipsoid, sum_squares, trid],
    "steep_ridges_and_drops": [dejong5, easom, michalewicz],
    "plate_shaped": [booth, matyas, mccormick, power_sum, zakharov],
    "valley_shaped": [camel3, camel6, dixon_price, rosenbrock],
    "other": [beale, branin, branin_modified, goldstein_price, permdb, styblinski_tang],
}
