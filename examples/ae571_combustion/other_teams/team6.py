"""AE 571 final project, team 6 (``Gasoline``, ``Hydrogen``, ``atombalance``, ``propertycalculator``;
O. Gonzalez and J. Oswald): dual (compression-ignition) cycle on gasoline or hydrogen.

Method: heat of combustion at T2 per kilogram of mixture; the heat input
as a function of beta,
``q(beta) = [H_R(T3) - H_R(T2)] + [H_P(beta T3) - H_P(T3)] - R_R (T3 - T2)``,
with the constant-volume step on the *reactant* enthalpy and the
constant-pressure step on the products, set equal to the heating value;
``v1`` from the gas constant of air. Plots of q(beta) against the heating
value and against the cycle work, as in the MATLAB.

Fixes: the beta search stepped by 1e-2 down to 1e-10 for up to 1000
iterations and rounded to two decimals (a root finder here); the
"required mass flow" ``55 / W`` is the mixture flow, the fuel flow is
``P / (eta LHV_fuel)``. The rich-mixture (CO) branch of the atom balance
was empty in the MATLAB; only lean mixtures are handled. The peak-pressure
prompt is ``--p-max``.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from unicodes.thermo import atom_balance, enthalpy_molar
from unicodes.thermo.engine_cycles import MOLAR_MASS, dual_cycle_work

RU = 8.3145


def setup(fuel, eps, er, T1, p1, n1, n2, p3):
    rxn = atom_balance(fuel, er)
    sp, Mf = rxn.fuel.species, rxn.fuel.molar_mass
    T2 = T1 * eps ** (n1 - 1)
    react = {sp: 1.0, "O2": rxn.a, "N2": rxn.n2}
    prod = {"CO2": rxn.b, "H2O": rxn.c, "O2": rxn.d, "N2": rxn.n2}

    def H(moles, T):
        return sum(n * enthalpy_molar(s, T) for s, n in moles.items())

    hc = H(react, T2) - H(prod, T2)  # J per mol fuel
    m = sum(n * MOLAR_MASS[s] for s, n in prod.items()) / 1000
    n_r = sum(react.values())
    R_R = RU * n_r / (Mf + rxn.a * 32 + rxn.n2 * 28) * 1000
    alpha = p3 / (p1 * eps**n1)
    T3 = alpha * T2

    def qin(beta):  # per kg of mixture
        return (H(react, T3) - H(react, T2) + H(prod, beta * T3) - H(prod, T3)) / m - R_R * (T3 - T2)

    return rxn, hc, m, alpha, qin


def run(fuel="gasoline (heavy)", eps=16.0, er=0.8, T1=298.0, p1=100e3, n1=1.38, n2=1.25, p3=5.5e6, power=55e3):
    rxn, hc, m, alpha, qin = setup(fuel, eps, er, T1, p1, n1, n2, p3)
    hc_mix = hc / m
    beta = brentq(lambda b: qin(b) - hc_mix, 1.0, 10.0)
    v1 = RU / 28.9647e-3 * T1 / p1
    w = dual_cycle_work(p1, v1, eps, alpha, beta, n1, n2)
    eta = w / hc_mix
    return dict(lhv=hc_mix, lhv_fuel=hc / (rxn.fuel.molar_mass / 1000), alpha=alpha, beta=beta, work=w,
                mep=w / (v1 - v1 / eps), eta=eta, fuel_flow=power / (eta * hc / (rxn.fuel.molar_mass / 1000)),
                v1=v1, qin=qin, rxn=rxn)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--fuel", choices=["gasoline (heavy)", "hydrogen"], default="gasoline (heavy)")
    p.add_argument("--p-max", type=float, default=5.5e6)
    a = p.parse_args()
    if a.p_max < 100e3 * 16**1.38:
        raise SystemExit(f"ERROR: P3 must be >= P2 = {100e3 * 16**1.38:.2f} Pa")
    r = run(a.fuel, p3=a.p_max)
    print("CYCLE DESIGN PARAMETERS: epsilon 16, ER 0.8, T1 298 K, P1 100 kPa, n1 1.38, n2 1.25")
    print(f"LOWER HEATING VALUES: fuel {r['lhv_fuel'] / 1e3:.2f} kJ/kg_fuel, mixture {r['lhv'] / 1e3:.2f} kJ/kg_mix")
    print(f"PV DIAGRAM VARIABLES: alpha {r['alpha']:.2f}, beta {r['beta']:.2f}")
    print(f"NET WORK {r['work'] / 1e3:.2f} kJ/kg, MEP {r['mep'] / 1e3:.2f} kPa, THERMAL EFFICIENCY {100 * r['eta']:.2f} %, "
          f"FUEL FLOW {r['fuel_flow']:.6f} kg/s")
    b = np.linspace(1.0, 5.0, 100)
    q = np.array([r["qin"](x) for x in b]) / 1e3
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    ax[0].plot(b, q, "b", label="qin")
    ax[0].axhline(r["lhv"] / 1e3, color="r", label="hc")
    ax[0].set(title="hC & qIN over beta", xlabel="beta", ylabel="Energy (kJ/kg)")
    ax[0].legend()
    ax[1].plot(b, q, "b", label="qin")
    ax[1].plot(b, [dual_cycle_work(100e3, r["v1"], 16, r["alpha"], x, 1.38, 1.25) / 1e3 for x in b], "r", label="Work")
    ax[1].set(title="Work & qIN over beta", xlabel="beta", ylabel="Energy (kJ/kg)")
    ax[1].legend()
    plt.show()
