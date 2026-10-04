"""AE 571 final project, group 3 (``Group3_AE571_Final_Project``): dual cycle for heavy gasoline,
light gasoline or hydrogen, with a P-v diagram.

Method: heat of combustion at T2 per kilogram of products. If releasing it
all at constant volume stays below the peak pressure, the cycle is an Otto
cycle (beta = 1); otherwise T3 = alpha T2 and the remaining heat raises the
products at constant pressure to T4 (temperatures found by iterating the
product enthalpies). ``v1`` from the gas constant of air.

Fixes: the fuel enthalpy used 4148 for 4184 J/kcal; T5 used ``n1``
instead of ``n2``; the fuel flow was ``55 / (eta W)`` (kg/s of nothing in
particular) and is ``P / (eta LHV_fuel)`` here. The temperature iteration
(steps of ``10 |error|`` K) is replaced by a root finder. The fuel and
peak-pressure prompts are options.
"""

import argparse

import matplotlib.pyplot as plt
from scipy.optimize import brentq

from unicodes.thermo import atom_balance, enthalpy_molar
from unicodes.thermo.engine_cycles import MOLAR_MASS, dual_cycle_pv, dual_cycle_work

RU = 8.314


def run(fuel="gasoline (heavy)", eps=16.0, er=0.8, T1=298.0, p1=100e3, n1=1.38, n2=1.25, p3=5.5e6, power=55e3):
    rxn = atom_balance(fuel, er)
    T2 = T1 * eps ** (n1 - 1)
    p2 = p1 * eps**n1
    prod = {"CO2": rxn.b, "H2O": rxn.c, "O2": rxn.d, "N2": rxn.n2}
    HR = enthalpy_molar(rxn.fuel.species, T2) + rxn.a * enthalpy_molar("O2", T2)  # N2 cancels
    HP = sum(n * enthalpy_molar(s, T2) for s, n in prod.items() if s != "N2")
    lhv_fuel = (HR - HP) / (rxn.fuel.molar_mass / 1000)
    n_mix = sum(prod.values())
    M_av = sum(n * MOLAR_MASS[s] for s, n in prod.items()) / n_mix
    lhv = (HR - HP) / (n_mix * M_av / 1000)

    def rise(T0, T, R):  # J/kg of products from T0 to T, minus R (T - T0) for constant volume
        dh = sum(n / n_mix * (enthalpy_molar(s, T) - enthalpy_molar(s, T0)) for s, n in prod.items())
        return (dh - R * (T - T0)) / (M_av / 1000)

    T_test = brentq(lambda T: rise(T2, T, RU) - lhv, T2, 6000)
    if T_test / T2 <= p3 / p2:
        alpha, beta, T3 = T_test / T2, 1.0, T_test
    else:
        alpha = p3 / p2
        T3 = alpha * T2
        T4 = brentq(lambda T: rise(T3, T, 0.0) - (lhv - rise(T2, T3, RU)), T3, 6000)
        beta = T4 / T3
    v1 = 8.313 / 28.96e-3 * T1 / p1
    w = dual_cycle_work(p1, v1, eps, alpha, beta, n1, n2)
    eta = w / lhv
    return dict(lhv=lhv, alpha=alpha, beta=beta, work=w, mep=w / (v1 - v1 / eps), eta=eta,
                fuel_flow=power / (eta * lhv_fuel), v1=v1, T3=T3, rxn=rxn)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--fuel", choices=["gasoline (heavy)", "gasoline (light)", "hydrogen"], default="gasoline (heavy)")
    p.add_argument("--p-max", type=float, default=5.5e6)
    a = p.parse_args()
    r = run(a.fuel, p3=a.p_max)
    print("fuel-lean mixture (ER < 1)")
    for k, unit, scale in (("lhv", "kJ/kg", 1e-3), ("alpha", "", 1), ("beta", "", 1), ("work", "kJ/kg", 1e-3),
                           ("mep", "kPa", 1e-3), ("eta", "", 1), ("fuel_flow", "kg/s", 1)):
        print(f"  {k:10s} {r[k] * scale:12.5g} {unit}")
    v, pv = dual_cycle_pv(100e3, r["v1"], 16, r["alpha"], r["beta"], 1.38, 1.25)
    plt.plot(v, pv / 1e3)
    plt.xlabel("Specific Volume, v [m^3/kg]")
    plt.ylabel("Pressure, P [kPa]")
    plt.title("P vs v")
    plt.grid(True)
    plt.show()
