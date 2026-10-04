"""AE 700 project 1: visible push-broom imager sizing (``project_1_v_1_0``, ``Project1Runthrough_5``,
``optimized_project_1_V_2``).

The interactive ``input()`` prompts become command-line arguments. The
design is checked against the course requirements (resolution <= 3.5 m,
aperture <= 0.1 m, swath >= 6 mi, 1 m <= f <= 10 m, signal >= 0.3 V/um,
aircraft below 20 km or satellite above 160 km). ``--optimize`` runs the
SLSQP version of ``optimized_project_1_V_2`` (``fmincon``) on
``x = [D, SW, R]``. That formulation was unfinished in the MATLAB (the
signal target and the signal expression are in different units, so the
optimum sits on the bounds); it is reproduced as written, with bounds added.
"""

import argparse

import numpy as np
from scipy.optimize import minimize

from unicodes import remote_sensing as rs
from unicodes.atmosphere import isa

p = argparse.ArgumentParser()
p.add_argument("--time", type=float, default=12.0, help="time of day, decimal 24 h (6-18)")
p.add_argument("--resolution", type=float, default=3.5)
p.add_argument("--aperture", type=float, default=0.1)
p.add_argument("--focal-length", type=float, default=2.0)
p.add_argument("--optimize", action="store_true")
a = p.parse_args()

lam = 0.5e-6
w_d = 3000 * 60e-6  # 3000 detectors of 60 um
A_det = 60e-6**2
responsivity = 1e10
reflectance, E_sun, tau = 0.3, 1370.0, 0.96

if not 6 <= a.time <= 18:
    raise SystemExit("Time of day is outside of time range")
theta = rs.sun_zenith_from_time(a.time)
if a.resolution > 3.5:
    raise SystemExit("Resolution is too coarse")
if not 0 < a.aperture <= 0.1:
    raise SystemExit("Aperture diameter out of range")
H = rs.altitude_for_resolution(a.resolution, lam, a.aperture)
kind = "satellite" if H >= 160e3 else "aircraft" if H <= 20e3 else None
print(f"Altitude {H:.0f} m:", kind or "not at a functional height")
if not 1 <= a.focal_length <= 10:
    raise SystemExit("Focal length out of range")
SW = rs.swath_width(H, w_d, a.focal_length)
print(f"Swath width {SW:.0f} m", "(too narrow)" if SW < 9656.064 else "")
L = rs.reflected_radiance(E_sun, reflectance, theta)
E = rs.image_irradiance(L, np.pi / 4 * a.aperture**2, a.focal_length, rs.rayleigh_angle(lam, a.aperture), tau)
S = rs.detector_signal(responsivity, E, A_det)
print(f"Reflected radiance {L:.2f} W/m^2/sr, sensor irradiance {E:.4g} W/m^2, signal {S:.4g}",
      "(too low)" if S < 0.3 else "")

if a.optimize:
    S_min = 300.0  # mV/um as in the MATLAB objective
    theta_dot_sun = np.deg2rad(15 / 3600)
    A_search = 23567e3
    mu_earth, R_earth = 6.674e-11 * 5.972e24, 6371e3

    def speed(H):
        if H <= 20e3:
            return 0.8 * isa(H).a
        if H >= 160e3:
            return np.sqrt(mu_earth / (R_earth + H))
        return np.nan  # no vehicle between 20 and 160 km (the MATLAB gave V = 0)

    def objective(x):
        D, SW, R = x
        H = rs.altitude_for_resolution(R, lam, D)
        V = speed(H)
        sun = np.cos(theta_dot_sun * A_search / (2 * SW * V))
        signal = 1.22**2 * tau * w_d**2 * lam * np.cos(1.22 * lam / D) ** 4 * E_sun * reflectance * sun / (4 * SW**2 * R**2)
        out = (signal - S_min / (responsivity * A_det)) ** 2
        return out if np.isfinite(out) else 1e30

    def constraints(x):  # >= 0
        D, SW, R = x
        f = SW * R * D / (1.22 * w_d * lam)
        return np.array([D - 0.1e-3, 0.1 - D, SW - 6 * 1609.344, 3.5 - R, f - 1, 10 - f])

    # Bounds added: without them the MATLAB formulation drives R (and H) to zero.
    res = minimize(objective, [0.05, 6 * 1609.344, 3.5], method="SLSQP", bounds=[(1e-4, 0.1), (9656.064, 1e6), (0.5, 3.5)],
                   constraints={"type": "ineq", "fun": constraints})
    D, SW, R = res.x
    print(f"Optimised: D = {D:.4g} m, SW = {SW:.0f} m, R = {R:.3g} m, "
          f"H = {rs.altitude_for_resolution(R, lam, D):.0f} m, objective {res.fun:.3g}")
