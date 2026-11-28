"""AE 571 final project, group 9 (``AE571_FinalProject_Group9``): dual cycle on gasoline or hydrogen.

Method: heat of combustion at T2 per kilogram of products; T3 = alpha T2;
a mean product cp (and cv = cp - R) from the enthalpy change between T2
and T3, used for both the constant-volume and the constant-pressure heat
release; ``v1`` from the gas constant of air.

Fixes: the gasoline enthalpy coefficient a1 was typed -24.978 (Heywood:
-24.078); "fuel required" divided the power by the cycle work, which gives
the mixture flow; the fuel flow is ``P / (eta LHV_fuel)``. The prompts for
fuel and manual parameters are options.
"""

import argparse

from unicodes.thermo import atom_balance, enthalpy_molar
from unicodes.thermo.engine_cycles import MOLAR_MASS, dual_cycle_work

RU = 8.314


def run(fuel="gasoline (heavy)", eps=16.0, er=0.8, T1=298.0, p1=100e3, n1=1.38, n2=1.25, p3=5.5e6, power=55e3):
    rxn = atom_balance(fuel, er)
    sp = rxn.fuel.species
    T2 = T1 * eps ** (n1 - 1)
    prod = {"CO2": rxn.b, "H2O": rxn.c, "O2": rxn.d, "N2": rxn.n2}
    HR = enthalpy_molar(sp, T2) + rxn.a * enthalpy_molar("O2", T2) + rxn.n2 * enthalpy_molar("N2", T2)
    HP = sum(n * enthalpy_molar(s, T2) for s, n in prod.items())
    m_prod = sum(n * MOLAR_MASS[s] for s, n in prod.items()) / 1000  # kg per mol fuel
    lhv = (HR - HP) / m_prod  # J/kg products
    p2 = p1 * eps**n1
    alpha = p3 / p2
    T3 = alpha * T2
    cp_i = {s: (enthalpy_molar(s, T3) - enthalpy_molar(s, T2)) / (T3 - T2) for s in prod}
    cp = sum(n * cp_i[s] for s, n in prod.items()) / m_prod
    cv = sum(n * (cp_i[s] - RU) for s, n in prod.items()) / m_prod
    q1 = cv * (T3 - T2)
    T4 = T3 + (lhv - q1) / cp
    beta = T4 / T3
    v1 = RU / 0.029 * T1 / p1
    w = dual_cycle_work(p1, v1, eps, alpha, beta, n1, n2)
    eta = w / lhv
    return dict(lhv=lhv, alpha=alpha, beta=beta, work=w, mep=w / (v1 - v1 / eps), eta=eta,
                fuel_flow=power / (eta * (HR - HP) / (rxn.fuel.molar_mass / 1000)), rxn=rxn)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--fuel", choices=["gasoline (heavy)", "hydrogen"], default="gasoline (heavy)")
    p.add_argument("--p-max", type=float, default=5.5e6)
    a = p.parse_args()
    r = run(a.fuel, p3=a.p_max)
    x = r["rxn"]
    print(f"{x.fuel.name} + {x.a:.4f} (O2 + 3.76 N2) -> {x.b:g} CO2 + {x.c:g} H2O + {x.d:.4f} O2 + {x.n2:.4f} N2")
    print(f"LHV {r['lhv'] / 1e3:.2f} kJ/kg, alpha {r['alpha']:.4f}, beta {r['beta']:.4f}, work {r['work'] / 1e3:.2f} kJ/kg, "
          f"mep {r['mep'] / 1e3:.1f} kPa, efficiency {r['eta']:.4f}, fuel {r['fuel_flow'] * 1e3:.4f} g/s")
