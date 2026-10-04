"""Elevator doublet response of the Cessna 182, rigid and with passive aero-compliant flaps.

Covers ``BenMays6DOF``, ``BenMays6DOF_PAH``, ``BenMays6DOF_PAH_sailplane(_1Flap/_2Flap)``,
``BenMays6DOF_PAH_cessna182_2Flap``, ``Cessna_182_PAH`` and ``Roskam_6DOF_C182``:
the same 6-DOF model with zero, one or two PAH flap degrees of freedom.
"""

import matplotlib.pyplot as plt
import numpy as np

from unicodes.flight_dynamics import sixdof as s6

ac = s6.cessna_182()
dt, T, angle = 0.01, 10.0, np.deg2rad(5)
t = np.arange(0, T + dt / 2, dt)
controls = np.zeros((len(t), 4))
controls[int(1 / dt):int(1.5 / dt), 1] = angle  # elevator pulse

flaps = s6.wing_flaps(ac.b / 2, 0.25, tau=0.45, I_h=0.9, k=0.2, c_r=0.05, delta_0=np.deg2rad(5))
runs = {
    "rigid": s6.simulate(ac, t, s6.initial_state(ac, theta=0.0), controls),
    "PAH": s6.simulate(ac, t, s6.initial_state(ac, theta=0.0, flaps=flaps), controls, flaps),
}
fig, axes = plt.subplots(4, 1, sharex=True, figsize=(7, 9))
for name, r in runs.items():
    axes[0].plot(t, np.rad2deg(r["alpha"]), label=name)
    axes[1].plot(t, np.rad2deg(r["Q"]))
    axes[2].plot(t, r["z_D"])
for name in ("delta_f1", "delta_f2"):
    axes[3].plot(t, np.rad2deg(runs["PAH"][name]), label=name)
for ax, label in zip(axes, [r"$\alpha$ (deg)", "Q (deg/s)", "D (ft)", r"$\delta_f$ (deg)"]):
    ax.set_ylabel(label)
    ax.grid(True)
axes[0].legend()
axes[3].legend()
axes[-1].set_xlabel("time (s)")
plt.show()
