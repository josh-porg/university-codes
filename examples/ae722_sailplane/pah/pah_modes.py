"""Longitudinal eigenmodes of the Cessna 182 with a PAH flap (``BenMays_6DOF_jacobian``).

Linearises the full model by central differences, removes position,
actuator, lateral-directional and second-flap states, and prints the
eigenvalues with the dominant states of each mode.
"""

import numpy as np

from unicodes.flight_dynamics import sixdof as s6

ac = s6.cessna_182()
flaps = s6.wing_flaps(ac.b / 2, 0.25, tau=0.45, I_h=0.01, k=0.2, c_r=0.05, delta_0=0.0)
x0 = s6.initial_state(ac, theta=0.0, flaps=flaps)
J = s6.linearize(ac, x0, np.zeros(4), flaps)
names = s6.STATE_NAMES + ["delta_f1", "delta_f1_dot", "delta_f2", "delta_f2_dot"]
drop = set(range(9, 16)) | {2, 3, 5, 6, 8, 18, 19}
keep = [i for i in range(len(names)) if i not in drop]
J = J[np.ix_(keep, keep)]
lam, W = np.linalg.eig(J)
for i in np.argsort(-lam.real):
    w = np.abs(W[:, i])
    top = [names[keep[j]] for j in np.argsort(-w)[:3]]
    print(f"lambda = {lam[i]:.4f}   dominant states: {', '.join(top)}")
