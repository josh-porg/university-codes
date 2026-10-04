"""Propeller blade optimisation and minimum-power radius study (AE 722 Propeller Optimization).

``CE_PrOptimization`` + ``PropDesign_SetPower``: for each flight speed and
radius, the cross-entropy method tunes the twist law and root chord to
maximise efficiency at 19.7 kW, rejecting designs whose outer sections stall.
``Plot_PrOptiamal_Solution``: fits eta_p(V) per radius and finds the
speed that minimises climb power at 800 fpm. Sample counts are reduced.
"""

import numpy as np
from scipy.optimize import curve_fit, minimize_scalar

from unicodes.aero import PropellerBlade, performance as perf
from unicodes.atmosphere import RHO0, isa
from unicodes.optimize import cross_entropy_maximize
from unicodes.units import FPM, FT, HP, KT, LBF, PSF

rho = float(isa(3000 * FT, 10).rho)
P_DESIRED = 19.7e3


def best_efficiency(V, radius, rng=0):
    def eta(x):
        blade = PropellerBlade(*x, radius=radius)
        try:
            p = blade.match_power(P_DESIRED, V=V, rho=rho)
        except ValueError:
            return 0.0
        return 0.0 if np.any(np.rad2deg(p.alpha[39:]) > 16.9) else p.efficiency

    x0 = [np.deg2rad(45), np.deg2rad(15), np.deg2rad(-42), -0.7, 0.1]
    res = cross_entropy_maximize(eta, x0, [0.06, 0.03, 0.03, 0.03, 0.06], n_samples=20, max_iterations=4, rng=rng)
    return res.fun


speeds = np.arange(40, 81, 10) * KT
radii = np.array([1.0, 1.2])
eta = np.array([[best_efficiency(V, r) for r in radii] for V in speeds])
print("eta_p [V x radius]:\n", np.round(eta, 3))

# Minimum climb power vs radius
A, e, C_D0, W, ws = 29.0, 0.9, 0.009, 700 * 9.81, 500.0
sigma = rho / RHO0
for j, r in enumerate(radii):
    (a, b), _ = curve_fit(lambda x, a, b: (x - a) / (x - b), speeds, eta[:, j], p0=(1.0, 0.0), maxfev=20000)
    eta_fit = lambda V: (V - a) / (V - b)  # noqa: E731

    def power(V):
        C_L = 2 * ws / (rho * V**2)
        index = C_L**1.5 / (C_D0 + C_L**2 / (np.pi * A * e))
        return W / LBF / eta_fit(V) * (800 / 33000 + np.sqrt(ws / PSF) / (19 * index * np.sqrt(sigma))) * HP

    res = minimize_scalar(power, bounds=(20, 50), method="bounded")
    print(f"radius {r:.1f} m: minimum climb power {res.fun / 1e3:.1f} kW at {res.x:.1f} m/s")
print("(check) Roskam climb power at 30 m/s, eta 0.85:",
      perf.power_for_climb(800 * FPM, 0.85, W, ws, sigma, C_D0, A, e) / 1e3, "kW")
