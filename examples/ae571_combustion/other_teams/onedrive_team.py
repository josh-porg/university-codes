"""AE 571 final project, team in ``OneDrive_1_12-16-2022`` (``AE571_FinalProject_PartI``/``PartII``,
``atombalancefinal``, ``enthalpycalculator2``, ``enthalpycalculatorfinal``, their ``H2`` variants and
``specificheats``): dual cycle, part I heavy gasoline, part II hydrogen.

Method: heat of combustion per kilogram of products; T3 = alpha T2;
``beta = 1 + (q - cv (alpha - 1) T2) / (cp alpha T2)`` with mean cp and cv
of the fuel-air (reactant) mixture; ``v1`` with R = 0.2891 kJ/(kg K).

Fixes: the reactant enthalpies were taken at T1 and the products at T2, so
the "heating value" included the compression heating (both at T2 here);
the hydrogen atom balance used y = 1 (H instead of H2); the hard-coded
"cp" and "cv" (1.1e4 to 1.6e5) were molar enthalpies printed by
``specificheats``, not heat capacities, so the reactant cp and cv here come
from the species polynomials averaged over the heat-release temperatures;
the work had ``CR^n1 - 1`` for ``CR^(n1-1)``; the fuel rate
``eta * 55 / 100 / q`` is ``P / (eta LHV_fuel)``.
"""

import numpy as np

from unicodes.thermo import atom_balance, cp_molar, enthalpy_molar
from unicodes.thermo.engine_cycles import MOLAR_MASS, dual_cycle_work

RU = 8.314


def run(fuel="gasoline (heavy)", eps=16.0, er=0.8, T1=298.0, p1=100e3, n1=1.38, n2=1.25, p3=5.5e6, power=55e3):
    rxn = atom_balance(fuel, er)
    sp, Mf = rxn.fuel.species, rxn.fuel.molar_mass
    T2 = T1 * eps ** (n1 - 1)
    react = {sp: 1.0, "O2": rxn.a, "N2": rxn.n2}
    prod = {"CO2": rxn.b, "H2O": rxn.c, "O2": rxn.d, "N2": rxn.n2}
    dH = sum(n * enthalpy_molar(s, T2) for s, n in react.items()) - sum(n * enthalpy_molar(s, T2) for s, n in prod.items())
    m = sum(n * MOLAR_MASS[s] for s, n in prod.items()) / 1000
    q = dH / m
    alpha = p3 / (p1 * eps**n1)
    T3 = alpha * T2
    m_r = (Mf + rxn.a * MOLAR_MASS["O2"] + rxn.n2 * MOLAR_MASS["N2"]) / 1000
    Ts = np.linspace(T2, 2.5 * T3, 50)  # mean over the heat-release range
    cp = np.mean([sum(n * cp_molar(s, T) for s, n in react.items()) for T in Ts]) / m_r
    cv = cp - RU * sum(react.values()) / m_r
    beta = 1 + (q - cv * (alpha - 1) * T2) / (cp * alpha * T2)
    v1 = 289.1 * T1 / p1
    w = dual_cycle_work(p1, v1, eps, alpha, beta, n1, n2)
    eta = w / q
    return dict(lhv=q, alpha=alpha, beta=beta, work=w, mep=w / (v1 - v1 / eps), eta=eta,
                fuel_flow=power / (eta * dH / (Mf / 1000)), rxn=rxn)


if __name__ == "__main__":
    for part, fuel in (("I", "gasoline (heavy)"), ("II", "hydrogen")):
        r = run(fuel)
        print(f"Part {part} ({fuel}): delta_hc {r['lhv'] / 1e3:.2f} kJ/kg, alpha {r['alpha']:.4f}, beta {r['beta']:.4f}, "
              f"W_cycle {r['work'] / 1e3:.2f} kJ/kg, mep {r['mep'] / 1e3:.1f} kPa, eta_th {r['eta']:.4f}, "
              f"fuel rate {r['fuel_flow'] * 1e3:.4f} g/s")
