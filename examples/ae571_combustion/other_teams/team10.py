"""AE 571 final project, team 10 (``Part_I``, ``Part_II``): dual cycle on gasoline (I) and hydrogen (II).

Method: heat of combustion at T2 per kilogram of products; T3 = alpha T2;
beta is raised in steps of 0.001 until the heat taken up by the fuel-air
(reactant) mixture, ``int cv dT`` from T2 to T3 plus ``int cp dT`` from
T3 to T4, reaches the heating value; ``v1`` from the gas constant of air.

Fixes: the reactant cv was ``int cp dT - R`` instead of
``int cp dT - R (T3 - T2)``; the fuel cp was integrated over theta = T/1000,
which divides its heat by 1000; the O2/N2 integrals were split at 1000 K
assuming T2 < 1000 K < T3 (the polynomials choose the range here); the
fuel flow ``55 / LHV_fuel`` is ``P / (eta LHV_fuel)``. The beta search is a
root finder. The peak-pressure prompt is ``--p-max``.
"""

import argparse

from scipy.integrate import quad
from scipy.optimize import brentq

from unicodes.thermo import atom_balance, cp_molar, enthalpy_molar
from unicodes.thermo.engine_cycles import MOLAR_MASS, dual_cycle_work

RU = 8.3145


def run(fuel="gasoline (heavy)", eps=16.0, er=0.8, T1=298.0, p1=100e3, n1=1.38, n2=1.25, p3=5.5e6, power=55e3):
    rxn = atom_balance(fuel, er)
    sp, Mf = rxn.fuel.species, rxn.fuel.molar_mass
    T2 = T1 * eps ** (n1 - 1)
    react = {sp: 1.0, "O2": rxn.a, "N2": rxn.n2}
    prod = {"CO2": rxn.b, "H2O": rxn.c, "O2": rxn.d, "N2": rxn.n2}
    HC = sum(n * enthalpy_molar(s, T2) for s, n in react.items()) - sum(n * enthalpy_molar(s, T2) for s, n in prod.items())
    hcf = HC / (Mf / 1000)
    hc = HC / (sum(n * MOLAR_MASS[s] for s, n in prod.items()) / 1000)
    alpha = p3 / (p1 * eps**n1)
    T3 = alpha * T2
    m_r = (Mf + rxn.a * 32 + rxn.n2 * 28) / 1000
    n_r = sum(react.values())

    def cp_r(T):
        return sum(n * cp_molar(s, T) for s, n in react.items())

    q23 = (quad(cp_r, T2, T3)[0] - n_r * RU * (T3 - T2)) / m_r

    def heat(beta):
        return q23 + quad(cp_r, T3, beta * T3)[0] / m_r - hc

    beta = brentq(heat, 1.0, 10.0)
    v1 = RU / 28.9647e-3 * T1 / p1
    w = dual_cycle_work(p1, v1, eps, alpha, beta, n1, n2)
    eta = w / hc
    return dict(lhv=hc, alpha=alpha, beta=beta, work=w, mep=w / (v1 - v1 / eps), eta=eta,
                fuel_flow=power / (eta * hcf), rxn=rxn)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--p-max", type=float, default=5.5e6)
    a = p.parse_args()
    for part, fuel in (("I", "gasoline (heavy)"), ("II", "hydrogen")):
        r = run(fuel, p3=a.p_max)
        x = r["rxn"]
        print(f"Part {part}: 1) {x.fuel.name} + {x.a:.5f}(O2 + 3.76N2) = {x.b:.2f}CO2 + {x.c:.2f}H2O + {x.d:.3f}O2 + {x.n2:.3f}N2")
        print(f"  2) LHV at T2 {r['lhv'] / 1e3:.2f} kJ/kg  3) alpha {r['alpha']:.3f}, beta {r['beta']:.3f}  "
              f"4) work {r['work'] / 1e3:.2f} kJ/kg  5) mep {r['mep'] / 1e3:.2f} kPa  6) efficiency {r['eta']:.4f}  "
              f"7) {r['fuel_flow']:.6f} kg fuel/s")
