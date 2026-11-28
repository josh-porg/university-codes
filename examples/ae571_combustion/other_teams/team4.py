"""AE 571 final project, team 4 (``AE_571_Final_Project_Code``): dual cycle for methane, propane,
gasoline, hydrogen or diesel, or a custom CxHy fuel with Heywood coefficients.

Method: heat of combustion at the reference temperature (298 K) per
kilogram of products; T3 = alpha T2 with ``q23 = int cv_p dT`` (products),
the rest at constant pressure to T4 (``int cp_p dT``). Then
``v4 = R_p T4 / p3`` with the products' gas constant while ``v3 = v1 / eps``
uses R = 287, so ``beta = v4 / v3`` is not T4/T3; that choice is kept. When
v4 <= v3 the whole heat goes in at constant volume (beta = 1).

Fixes: the product O2 enthalpy was taken at T2 while everything else was
at 298 K; the fuel enthalpy was evaluated at T2 (theta = T2/1000) in the
same balance; the molar mass of water was 20.16; the fuel flow
``55 kW / LHV_fuel`` is ``P / (eta LHV_fuel)``. Prompts are options.
"""

import argparse

from scipy.integrate import quad
from scipy.optimize import brentq

from unicodes.thermo import FUELS, Fuel, atom_balance, cp_molar, enthalpy_molar
from unicodes.thermo.combustion import HEYWOOD_FUELS
from unicodes.thermo.engine_cycles import MOLAR_MASS, dual_cycle_work

RU = 8.31451
T_REF = 298.0


def run(fuel="gasoline (heavy)", eps=16.0, er=0.8, T1=T_REF, p1=100e3, n1=1.38, n2=1.25, p3=5.5e6, power=55e3):
    rxn = atom_balance(fuel, er)
    sp, Mf = rxn.fuel.species, rxn.fuel.molar_mass
    T2 = T1 * eps ** (n1 - 1)
    p2 = p1 * eps**n1
    prod = {"CO2": rxn.b, "H2O": rxn.c, "O2": rxn.d, "N2": rxn.n2}
    dH = (enthalpy_molar(sp, T_REF) + rxn.a * enthalpy_molar("O2", T_REF)
          - sum(n * enthalpy_molar(s, T_REF) for s, n in prod.items() if s != "N2"))
    delta_hc = dH / (Mf / 1000)
    m_p = sum(n * MOLAR_MASS[s] for s, n in prod.items()) / 1000
    lhv = dH / m_p
    R_mix = RU * sum(prod.values()) / m_p
    alpha = p3 / p2
    T3 = alpha * T2

    def cp(T):
        return sum(n * cp_molar(s, T) for s, n in prod.items()) / m_p

    v1 = 287 * T1 / p1
    v3 = v1 / eps
    q23 = quad(lambda T: cp(T) - R_mix, T2, T3)[0]
    T4 = brentq(lambda T: quad(cp, T3, T)[0] - (lhv - q23), T3, 6000)
    v4 = R_mix * T4 / p3
    if v4 <= v3:  # all heat at constant volume
        T3 = brentq(lambda T: quad(lambda t: cp(t) - R_mix, T2, T)[0] - lhv, T2, 6000)
        alpha, beta = T3 / T2, 1.0
    else:
        beta = v4 / v3
    w = dual_cycle_work(p1, v1, eps, alpha, beta, n1, n2)
    eta = w / lhv
    return dict(lhv=lhv, alpha=alpha, beta=beta, work=w, mep=w / (v1 - v3), eta=eta,
                fuel_flow=power / (eta * delta_hc), rxn=rxn)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--fuel", choices=["methane", "propane", "gasoline (heavy)", "hydrogen", "diesel"],
                   default="gasoline (heavy)")
    p.add_argument("--custom", type=float, nargs=9, metavar=("X", "Y", "MW", "A1", "A2", "A3", "A4", "A5", "A6"),
                   help="custom CxHy fuel with Heywood enthalpy coefficients")
    p.add_argument("--p-max", type=float, default=5.5e6)
    a = p.parse_args()
    fuel = a.fuel
    if a.custom:
        x, y, mw, *coef = a.custom
        key = f"C{x:g}H{y:g}"
        HEYWOOD_FUELS[key] = tuple(coef)
        fuel = Fuel("Custom", x, y, mw, key)
    r = run(fuel if isinstance(fuel, Fuel) else FUELS[fuel], p3=a.p_max)
    x = r["rxn"]
    print(f"C_{x.fuel.x:g}H_{x.fuel.y:g} + {x.a:.2f}(O2 + 3.76N2) -> {x.b:.2f}CO2 + {x.c:.2f}H2O + {x.d:.2f}O2 + {x.n2:.2f}N2")
    print(f"LHV {r['lhv'] / 1e3:.2f} kJ/kg, alpha {r['alpha']:.4f}, beta {r['beta']:.4f}, W_cycle {r['work'] / 1e3:.2f} kJ/kg, "
          f"mep {r['mep'] / 1e6:.4f} MPa, efficiency {r['eta']:.4f}, fuel {r['fuel_flow']:.6f} kg/s")
