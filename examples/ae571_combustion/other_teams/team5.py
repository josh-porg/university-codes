"""AE 571 final project, team 5 ("Corrected Code that works": ``AE571_FinalProject_Part1``,
``AE571_Part1_Calculators``, ``AE571_FinalProject_Part2_Complete``, ``AtomBalanceLab1``,
``enthalpycalculatorAir``/``Fuel``/``O2``/``Products``): dual cycle, part I gasoline, part II hydrogen.

Both parts compute the heat of combustion at T2 per kilogram of mixture
and the constant-pressure adiabatic flame temperature T_ad of the
reactants at 298 K, but close the cycle differently:

* Part I: ``beta = (q - (cv3 T3 - cv2 T2) + cp3 T3) / (cp4 T3)`` with
  product cp and cv evaluated at T2, T3 and T4 = T_ad;
* Part II: ``beta = v4 / v2`` with ``v4 = R_p T_ad / p3``; when the
  products at T_ad and v2 stay below the peak pressure, an Otto cycle.

``v1`` from the reactant gas constant.

Fixes (part I): the H2O, O2 and N2 product enthalpies had the bracket
closed after the first coefficient (``Ru*(a1)*T + ...``), so only one term
was multiplied by R; the fuel enthalpy used ``a5/5 theta^5`` for
``-a5/theta``; the flame-temperature scan in steps of 0.75 K stopped at
1 % (a root finder here); the "fuel required" divided by the mixture
heating value and is ``P / (eta LHV_fuel)`` in both parts.
"""

from unicodes.thermo import adiabatic_flame_temperature, atom_balance, cp_molar, enthalpy_molar
from unicodes.thermo.engine_cycles import MOLAR_MASS, dual_cycle_work

RU = 8.3145


def _common(fuel, eps, er, T1, p1, n1):
    rxn = atom_balance(fuel, er)
    sp, Mf = rxn.fuel.species, rxn.fuel.molar_mass
    T2 = T1 * eps ** (n1 - 1)
    react = {sp: 1.0, "O2": rxn.a, "N2": rxn.n2}
    prod = {"CO2": rxn.b, "H2O": rxn.c, "O2": rxn.d, "N2": rxn.n2}
    dH = sum(n * enthalpy_molar(s, T2) for s, n in react.items()) - sum(n * enthalpy_molar(s, T2) for s, n in prod.items())
    m = (Mf + rxn.a * MOLAR_MASS["O2"] + rxn.n2 * MOLAR_MASS["N2"]) / 1000
    v1 = RU * sum(react.values()) / m * T1 / p1
    return rxn, prod, T2, dH, m, v1, dH / (Mf / 1000)


def run_part1(fuel="gasoline (heavy)", eps=16.0, er=0.8, T1=298.0, p1=100e3, n1=1.38, n2=1.25, p3=5.5e6, power=55e3):
    rxn, prod, T2, dH, m, v1, lhv_fuel = _common(fuel, eps, er, T1, p1, n1)
    lhv = dH / m
    alpha = p3 / (p1 * eps**n1)
    T3 = alpha * T2
    T4 = adiabatic_flame_temperature(rxn.fuel, er, T_reactants=298.0)

    def cp(T):
        return sum(n * cp_molar(s, T) for s, n in prod.items()) / m

    def cv(T):
        return sum(n * (cp_molar(s, T) - RU) for s, n in prod.items()) / m

    beta = (lhv - (cv(T3) * T3 - cv(T2) * T2) + cp(T3) * T3) / (T3 * cp(T4))
    w = dual_cycle_work(p1, v1, eps, alpha, beta, n1, n2)
    eta = w / lhv
    return dict(lhv=lhv, alpha=alpha, beta=beta, work=w, mep=w / (v1 - v1 / eps), eta=eta,
                fuel_flow=power / (eta * lhv_fuel), T_ad=T4, rxn=rxn)


def run_part2(fuel="hydrogen", eps=16.0, er=0.8, T1=298.0, p1=100e3, n1=1.38, n2=1.25, p3=5.5e6, power=55e3):
    rxn, prod, T2, dH, m, v1, lhv_fuel = _common(fuel, eps, er, T1, p1, n1)
    lhv = dH / m
    T_ad = adiabatic_flame_temperature(rxn.fuel, er, T_reactants=298.0)
    v2 = v1 / eps
    R_p = RU * sum(prod.values()) / m
    p2 = p1 * eps**n1
    p_chem = R_p * T_ad / v2
    if p_chem > p3:
        alpha, beta = p3 / p2, R_p * T_ad / p3 / v2
    else:
        alpha, beta = p_chem / p2, 1.0
    w = dual_cycle_work(p1, v1, eps, alpha, beta, n1, n2)
    eta = w / lhv
    return dict(lhv=lhv, alpha=alpha, beta=beta, work=w, mep=w / (v1 - v2), eta=eta,
                fuel_flow=power / (eta * lhv_fuel), T_ad=T_ad, rxn=rxn)


run = run_part1

if __name__ == "__main__":
    for part, f, fuel in (("I", run_part1, "gasoline (heavy)"), ("II", run_part2, "hydrogen")):
        r = f(fuel)
        print(f"Part {part} ({fuel}): T_ad {r['T_ad']:.1f} K, LHV {r['lhv'] / 1e3:.2f} kJ/kg, alpha {r['alpha']:.4f}, "
              f"beta {r['beta']:.4f}, work {r['work'] / 1e3:.2f} kJ/kg, MEP {r['mep'] / 1e3:.1f} kPa, "
              f"efficiency {r['eta']:.4f}, fuel {r['fuel_flow']:.6f} kg/s")
