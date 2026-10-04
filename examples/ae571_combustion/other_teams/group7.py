"""AE 571 final project, group 7 (``AE_571_Group7_Final_Code``; M. Bonham, E. Horst, A. Kaushik,
B. Svoboda, J. Stockley, J. Wegiel): dual cycle, part I gasoline and part II hydrogen.

Method: heat of combustion at T2 per kilogram of products; ``v1`` from
R = 287; the heat input balance ``q = cp (T3 - T2) + cp (T4 - T3) - R (T3 - T2)``
solved for beta, with T3 = alpha T2.

Fixes: the "product cp" was ``|HP / M_p| / M_p`` (an enthalpy divided by a
molar mass twice); the balance is written here with the product enthalpies,
``q = h_p(T4) - h_p(T2) - R_p (T3 - T2)``. The gasoline enthalpy at T1
dropped the theta factors and the hydrogen enthalpy dropped its a2 term
(both replaced by the species polynomials); the work was scaled by 0.001
twice over (kJ and J mixed); the fuel flow assumed one mole of fuel per
cycle and the cycle work per mole, and is ``P / (eta LHV_fuel)`` here. The
maximum-pressure prompt is ``--p-max``.
"""

import argparse

from scipy.optimize import brentq

from unicodes.thermo import atom_balance, enthalpy_molar
from unicodes.thermo.engine_cycles import MOLAR_MASS, dual_cycle_work

RU = 8.3144598


def run(fuel="gasoline (heavy)", eps=16.0, er=0.8, T1=298.0, p1=100e3, n1=1.38, n2=1.25, p3=5.5e6, power=55e3):
    rxn = atom_balance(fuel, er)
    T2 = T1 * eps ** (n1 - 1)
    prod = {"CO2": rxn.b, "H2O": rxn.c, "O2": rxn.d, "N2": rxn.n2}
    HR = enthalpy_molar(rxn.fuel.species, T2) + rxn.a * enthalpy_molar("O2", T2) + rxn.n2 * enthalpy_molar("N2", T2)
    HP = sum(n * enthalpy_molar(s, T2) for s, n in prod.items())
    m_prod = sum(n * MOLAR_MASS[s] for s, n in prod.items()) / 1000
    lhv = (HR - HP) / m_prod
    R_mix = RU * sum(prod.values()) / m_prod
    p2 = p1 * eps**n1
    alpha = p3 / p2
    T3 = alpha * T2

    def h(T):
        return sum(n * enthalpy_molar(s, T) for s, n in prod.items()) / m_prod

    T4 = brentq(lambda T: h(T) - h(T2) - R_mix * (T3 - T2) - lhv, T3, 6000)
    beta = T4 / T3
    v1 = 287 * T1 / p1
    w = dual_cycle_work(p1, v1, eps, alpha, beta, n1, n2)
    eta = w / lhv
    return dict(lhv=lhv, alpha=alpha, beta=beta, work=w, mep=w / (v1 - v1 / eps), eta=eta,
                fuel_flow=power / (eta * (HR - HP) / (rxn.fuel.molar_mass / 1000)), rxn=rxn)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--p-max", type=float, default=5.5e6)
    a = p.parse_args()
    for part, fuel in (("I", "gasoline (heavy)"), ("II", "hydrogen")):
        r = run(fuel, p3=a.p_max)
        print(f"Part {part} ({fuel}): LHV {r['lhv'] / 1e3:.2f} kJ/kg, alpha {r['alpha']:.4f}, beta {r['beta']:.4f}, "
              f"W_net {r['work'] / 1e3:.2f} kJ/kg, mep {r['mep'] / 1e3:.1f} kPa, n_th {r['eta']:.4f}, "
              f"mdot {r['fuel_flow'] * 1e3:.4f} g/s")
