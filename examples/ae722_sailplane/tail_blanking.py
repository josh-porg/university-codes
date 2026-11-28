"""Horizontal-tail dynamic-pressure ratio from wind-tunnel blanking data (``Windtunnel_Blanking_Analysis_V0``,
``eta_h_Function_Example``).

Pass ``Blanking Data.mat`` (rows: alpha deg, -, V at the tail m/s); tunnel speed 9 m/s.
Saves ``Archytas_horizontal_Tail_Effectiveness.npz`` for the PAH moment study.
"""

import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

from unicodes.io import load_mat

d = np.asarray(load_mat(sys.argv[1])["BlankingData"], dtype=float)
alphas, V = d[0], d[2]
eta_h = V**2 / 9.0**2  # q_tail / q_freestream
np.savez("Archytas_horizontal_Tail_Effectiveness.npz", alphas=alphas, eta_h=eta_h)
fit = CubicSpline(alphas[:-1], eta_h[:-1])  # last point is bad
a = np.linspace(-15, 35, 200)
plt.scatter(alphas, eta_h, c=eta_h, cmap="cool")
plt.plot(a, fit(a), "g")
plt.xlabel(r"$\alpha$ (deg)")
plt.ylabel(r"$\eta_h$")
plt.show()
