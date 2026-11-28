"""Tune the PAH flap spring, dampers, neutral angle and inertia to minimise peak
acceleration of the AETHER jet in Dryden turbulence (``AETHER_Ben_6DOF_with_optimizer*``,
``Ben_6DOF_with_optimizer(_fast_script/_SA)``, ``Aether_6DOF``).

Objective (as in the MATLAB): the largest second difference of the
inertial position. The MATLAB GA picked parents with probability
proportional to the *objective*, which favours bad designs when minimising;
:func:`unicodes.optimize.genetic_minimize` uses the inverse. Population and
generation counts are small here; the MATLAB ran 50 x 20 in parallel.
"""

import numpy as np

from unicodes.flight_dynamics import sixdof as s6, turbulence as tb
from unicodes.optimize import genetic_minimize
from unicodes.units import FT

ac = s6.aether()
dt, T = 0.1, 10.0
t = np.arange(0, T + dt / 2, dt)
controls = np.zeros((len(t), 4))
controls[:, 1] = np.deg2rad(np.deg2rad(1.54635))  # trim elevator (the MATLAB applied deg2rad twice)
ug, vg, wg = tb.dryden_gusts(len(t), dt, ac.V_trim * FT, 46000 * FT, "severe", rng=0)
pg, qg, rg = tb.dryden_rotational_gusts(wg, vg, dt, ac.V_trim * FT, 46000 * FT, ac.b / 2 * FT, rng=1)
gusts = np.column_stack([ug / FT, 0 * vg, wg / FT, 0 * pg, qg, 0 * rg])  # lateral gusts off, as in V1_7

base = s6.wing_flaps(12.5416667, 4.42, tau=0.707184, I_h=0.9, k=0.2, c_r=0.05, delta_0=0.0, k_stop=1e6)


def simulate(x, disabled=False):
    """x = [k, c_r, c_c, delta_initial, I_h].

    As in the MATLAB, the fourth variable only sets the flaps' initial
    deflection (the spring neutral angle stayed 0) and Coulomb damping was
    forced to 0 inside the objective.
    """
    k, c_r, c_c, delta_init, I_h = x
    flaps = s6.with_pah(base, k=k, c_r=c_r, c_c=0.0, I_h=I_h, disabled=disabled)
    x0 = s6.initial_state(ac, flaps=flaps, elevator=controls[0, 1])
    x0[16::2] = delta_init
    return s6.simulate(ac, t, x0, controls, flaps, gusts, rtol=1e-5, atol=1e-7)


def objective(x):
    try:
        acc = simulate(x).inertial_acceleration()
    except Exception:
        return 1e9
    peak = np.max(np.linalg.norm(acc, axis=1))
    return peak if np.isfinite(peak) else 1e9


if __name__ == "__main__":
    rigid = simulate([0.2, 0.05, 0.0, 0.0, 0.9], disabled=True)
    print("rigid peak acceleration:", np.max(np.linalg.norm(rigid.inertial_acceleration(), axis=1)))
    lb = [-50, 0, 0.0004, -100, 1e-9]
    ub = [50, 50, 0.0005, 100, 50]
    res = genetic_minimize(objective, lb, ub, population_size=8, generations=3, n_elites=1,
                           mutation="uniform", mutation_rate=0.3, rng=0,
                           callback=lambda g, x, fx: print(f"generation {g}: best {fx:.4g}"))
    print("optimum [k, c_r, c_c, delta_0, I_h] =", res.x, "objective", res.fun)
