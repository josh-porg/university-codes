"""AE 571 final project, team 2 (``AE571_Group2_Final_Code``): dual cycle on gasoline and hydrogen.

Method: heat of combustion at T2 per kilogram of products; T3 = alpha T2;
T4 from the product enthalpy, ``h(T4) = h(T2) + q + R_p (T3 - T2)`` (the
constant-volume part costs ``R (T3 - T2)`` less than an enthalpy rise);
``v1`` from the reactant gas constant.

Fixes: the gasoline enthalpy added ``int cp dT`` from 298 K to an
enthalpy that already includes it; the T4 search stepped 0.01 K until the
enthalpies agreed within 5 kJ/kg (a root finder here); the fuel rate
``55 kW / LHV_fuel`` left out the efficiency and is ``P / (eta LHV_fuel)``.
"""

from scipy.optimize import brentq

from unicodes.thermo import atom_balance, enthalpy_molar
from unicodes.thermo.engine_cycles import dual_cycle_work

RU = 8.314
MW = {"C": 12, "H": 1, "O": 16, "N": 14}


def run(fuel="gasoline (heavy)", eps=16.0, er=0.8, T1=298.0, p1=100e3, n1=1.38, n2=1.25, p3=5.5e6, power=55e3):
    rxn = atom_balance(fuel, er)
    x, y, a = rxn.fuel.x, rxn.fuel.y, rxn.a
    Mf = x * MW["C"] + y * MW["H"]
    n_p = x + y / 2 + (a - x - y / 4) + 3.76 * a
    m_p = (x * 44 + y / 2 * 18 + (a - x - y / 4) * 32 + 3.76 * a * 28) / 1000  # kg per mol fuel
    M_r = (Mf + a * (32 + 3.76 * 28)) / (a * 4.76 + 1)
    prod = {"CO2": rxn.b, "H2O": rxn.c, "O2": rxn.d, "N2": rxn.n2}
    T2 = T1 * eps ** (n1 - 1)

    def Hp(T):
        return sum(n * enthalpy_molar(s, T) for s, n in prod.items())

    Hr = enthalpy_molar(rxn.fuel.species, T2) + a * (enthalpy_molar("O2", T2) + 3.76 * enthalpy_molar("N2", T2))
    dhc = (Hr - Hp(T2)) / (Mf / 1000)  # J/kg fuel
    lhv = (Hr - Hp(T2)) / m_p
    alpha = p3 / (p1 * eps**n1)
    T3 = alpha * T2
    R_p = RU * n_p / m_p
    target = lhv + Hp(T2) / m_p + R_p * (T3 - T2)
    T4 = brentq(lambda T: Hp(T) / m_p - target, T3, 6000)
    beta = T4 / T3
    v1 = RU / (M_r / 1000) * T1 / p1
    w = dual_cycle_work(p1, v1, eps, alpha, beta, n1, n2)
    eta = w / lhv
    return dict(lhv=lhv, alpha=alpha, beta=beta, work=w, mep=w / (v1 * (1 - 1 / eps)), eta=eta,
                fuel_flow=power / (eta * dhc), rxn=rxn)


if __name__ == "__main__":
    for fuel in ("gasoline (heavy)", "hydrogen"):
        r = run(fuel)
        print(f"{fuel}: LHV {r['lhv'] / 1e3:.2f} kJ/kg, alpha {r['alpha']:.3f}, beta {r['beta']:.3f}, "
              f"work {r['work'] / 1e3:.2f} kJ/kg, mep {r['mep'] / 1e3:.2f} kPa, efficiency {100 * r['eta']:.2f}%, "
              f"fuel {r['fuel_flow']:.6f} kg/s")
