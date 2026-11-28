"""Phase portrait and response of a single PAH flap (``adaptive_airfoil_dynamic_simulator``,
``hardening_spring_tester``, ``adaptive_airfoil_dynamics_linearizedEvaluation``)."""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from unicodes.flight_dynamics import sixdof as s6

flap = s6.PAHFlap(0.25, [0, -1, 0], [0, 1, 0], [-1, -1, 0], [0, 0, 0], I_h=1.0, k=2.0, c_r=10.0, c_c=2.0,
                  delta_0=np.deg2rad(10), C_h_alpha=-0.45, C_h_delta=-0.85, k_stop=1000.0)
rho, alpha, u = 1.225, np.deg2rad(2), 30.0
V = np.array([u, 0, u * np.tan(alpha)])
q = 0.5 * rho * V @ V
c_d = s6.critical_damping(flap.I_h, flap.k, q, flap.chord, flap.C_h_delta, zeta=0.7)
flap.c_r = max(c_d - flap.c_c, 0)
print(f"damping for zeta = 0.7: {c_d:.3f}; eigenvalues {s6.flap_eigenvalues(flap.I_h, flap.k, flap.c_r, q, flap.chord, flap.C_h_delta)}")


def rhs(t, y):
    return [y[1], flap.hinge_acceleration(y[0], y[1], rho, V, (0, 0, 0))]


sol = solve_ivp(rhs, (0, 5), [np.deg2rad(-10), 0], max_step=0.01)
d = np.deg2rad(20)
D, DD = np.meshgrid(np.linspace(-d, d, 15), np.linspace(-d, d, 15))
dD = DD
ddD = np.vectorize(lambda a, b: flap.hinge_acceleration(a, b, rho, V, (0, 0, 0)))(D, DD)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
ax1.contourf(np.rad2deg(D), np.rad2deg(DD), np.hypot(dD, ddD), 50, cmap="spring")
ax1.quiver(np.rad2deg(D), np.rad2deg(DD), dD, ddD)
ax1.plot(np.rad2deg(sol.y[0]), np.rad2deg(sol.y[1]), "k")
ax1.set(xlabel=r"$\delta$ (deg)", ylabel=r"$\dot\delta$ (deg/s)")
x = np.linspace(-2, 4, 400)
stop = s6.PAHFlap(1, [0, 0, 0], [0, 1, 0], [-1, 0, 0], [0, 0, 0], 1, 0, 0, 0, delta_min=-1, delta_max=3, k_stop=1, C_h_alpha=0, C_h_delta=0)
ax2.plot(x, [-stop.k_stop * (min(0, v - stop.delta_min) + max(0, v - stop.delta_max)) for v in x])
ax2.set(title="hardening stop spring", xlabel=r"$\delta$", ylabel="moment")
plt.show()
