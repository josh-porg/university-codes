"""Hohmann transfers, burn-to-orbit and rendezvous timing (``driver``, ``hoffman``,
``getTranserOrbitHoffmamn``, ``rendezvous``) and HW 5 (Kepler time of flight, J2 regression)."""

import numpy as np

from unicodes import orbital as o

RE = 6378e3
print("part 3: orbit after a prograde burn from 500 km")
for dv in (10, 100, 1000, 10000):
    ra, a, e = o.orbit_after_tangential_burn(6878e3, dv)
    print(f"  dv {dv:5d} m/s: r_a {ra / 1e3:12.1f} km, a {a / 1e3:12.1f} km, e {e:.4f}")
print("part 5:", o.hohmann(RE + 500e3, RE + 150e3))
print("part 7 (plane change only, 33 deg):", o.hohmann(RE + 130e3, RE + 130e3, np.deg2rad(90 - 57)))
print("part 8:", o.hohmann(RE + 150e3, RE + 20000e3, np.deg2rad(45 - 28)))
r = o.rendezvous(RE + 120e3, RE + 240e3, 0.0, np.deg2rad(135))
print(f"part 9: dv {r.transfer.dv_total:.2f} m/s, lead angle {np.rad2deg(r.lead_angle):.2f} deg, wait {r.wait_time / 3600:.2f} h")

nu = o.time_to_true_anomaly(0.0, 2.4e3, 6778e3, 0.04)
print(f"HW5: true anomaly after 2400 s = {np.rad2deg(nu):.3f} deg")
print(f"HW5: J2 nodal regression = {np.rad2deg(o.nodal_regression_rate(7e6, 0.1, np.pi / 4)) * 86400:.4f} deg/day")
