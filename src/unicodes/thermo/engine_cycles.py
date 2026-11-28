"""Fuel-air limited-pressure (dual) cycle with real-gas properties (AE 571 final project).

The cycle: polytropic compression 1-2 (``n1``), constant-volume heat release
2-3 up to the peak pressure, constant-pressure release 3-4, polytropic
expansion 4-5 (``n2``) and constant-volume rejection 5-1. With ``alpha =
p3 / p2``, ``beta = v4 / v3`` and ``delta = v5 / v4 = eps / beta`` the net
work per unit mass is (Ferguson)::

    w = p1 v1 eps^(n1-1) [alpha (beta - 1) + alpha beta / (n2 - 1) (1 - delta^(1-n2))
                          - (1 - eps^(1-n1)) / (n1 - 1)]

Heat is released into the products: ``u(T3) - u(T2)`` at constant volume,
then ``h(T4) - h(T3)`` at constant pressure, using the species
polynomials of :mod:`unicodes.thermo.combustion`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .combustion import R_UNIVERSAL, LeanCombustion, atom_balance, enthalpy_molar

MOLAR_MASS = {"CO2": 44.01, "H2O": 18.015, "O2": 31.999, "N2": 28.014, "CO": 28.010, "H2": 2.016}


def dual_cycle_work(p1, v1, eps, alpha, beta, n1, n2):
    """Net work per unit mass of the polytropic dual cycle (units of ``p1 v1``)."""
    delta = eps / beta
    return p1 * v1 * eps ** (n1 - 1) * (alpha * (beta - 1) + alpha * beta / (n2 - 1) * (1 - delta ** (1 - n2))
                                        - (1 - eps ** (1 - n1)) / (n1 - 1))


def dual_cycle_pv(p1, v1, eps, alpha, beta, n1, n2, n=100):
    """Closed ``(v, p)`` path of the cycle for a P-v diagram."""
    v2 = v1 / eps
    v4 = beta * v2
    vc = np.linspace(v1, v2, n)
    p2 = p1 * eps**n1
    p3 = alpha * p2
    ve = np.linspace(v4, v1, n)
    pe = p3 * (v4 / ve) ** n2
    v = np.concatenate([vc, [v2, v4], ve, [v1]])
    p = np.concatenate([p1 * (v1 / vc) ** n1, [p3, p3], pe, [p1]])
    return v, p


def products(rxn: LeanCombustion) -> dict:
    return {"CO2": rxn.b, "H2O": rxn.c, "O2": rxn.d, "N2": rxn.n2}


def reactants(rxn: LeanCombustion) -> dict:
    return {rxn.fuel.species: 1.0, "O2": rxn.a, "N2": rxn.n2}


def mixture_enthalpy(moles: dict, T) -> float:
    """``sum n_i h_i(T)`` in J (moles in mol)."""
    return sum(n * enthalpy_molar(s, T) for s, n in moles.items() if n)


def mixture_mass(moles: dict, fuel_molar_mass=None) -> float:
    """Mass in grams of a mixture given in moles (fuel species need ``fuel_molar_mass``)."""
    return sum(n * MOLAR_MASS.get(s, fuel_molar_mass or 0.0) for s, n in moles.items())


def heat_temperature(moles: dict, T0, heat, constant="v", T_max=6000.0) -> float:
    """Temperature reached when ``heat`` (J) goes into ``moles`` from ``T0`` at constant volume or pressure."""
    n_tot = sum(moles.values())
    h0 = mixture_enthalpy(moles, T0)

    def residual(T):
        dh = mixture_enthalpy(moles, T) - h0
        if constant == "v":
            dh -= n_tot * R_UNIVERSAL * (T - T0)
        return dh - heat

    return brentq(residual, T0, T_max)


@dataclass
class DualCycleResult:
    rxn: LeanCombustion
    lhv_fuel: float  # J/kg fuel
    lhv_mixture: float  # J/kg mixture
    alpha: float
    beta: float
    T: dict
    p: dict
    v1: float
    work: float  # J/kg mixture
    mep: float  # Pa
    efficiency: float
    fuel_flow: float  # kg/s for the requested power

    def summary(self) -> str:
        return (f"LHV {self.lhv_fuel / 1e6:.3f} MJ/kg fuel = {self.lhv_mixture / 1e6:.4f} MJ/kg mixture, "
                f"alpha = {self.alpha:.4f}, beta = {self.beta:.4f}\n"
                f"T2..T5 = {self.T[2]:.1f}, {self.T[3]:.1f}, {self.T[4]:.1f}, {self.T[5]:.1f} K; "
                f"net work {self.work / 1e3:.2f} kJ/kg, mep {self.mep / 1e6:.4f} MPa, efficiency {self.efficiency:.4f}, "
                f"fuel flow {self.fuel_flow * 1e3:.4f} g/s")


def dual_cycle(fuel="gasoline (heavy)", equivalence_ratio=0.8, eps=16.0, T1=298.0, p1=100e3, n1=1.38, n2=1.25,
               p_max=5.5e6, power=55e3) -> DualCycleResult:
    """The AE 571 project engine with the heat of combustion at T2 released into the products.

    If the whole heat at constant volume would stay below ``p_max``, the
    cycle is an Otto cycle (``beta = 1``) with the peak pressure that results.
    Mass basis: the fuel-air mixture (``v1`` from the reactant gas constant).
    """
    rxn = atom_balance(fuel, equivalence_ratio)
    R_mol = reactants(rxn)
    P_mol = products(rxn)
    Mf = rxn.fuel.molar_mass
    m_mix = mixture_mass(R_mol, Mf) / 1000  # kg per mol fuel
    T2 = T1 * eps ** (n1 - 1)
    p2 = p1 * eps**n1
    dH = mixture_enthalpy(R_mol, T2) - mixture_enthalpy(P_mol, T2)  # J per mol fuel at T2
    lhv_fuel = dH / (Mf / 1000)
    n_r = sum(R_mol.values())
    R_mix = n_r * R_UNIVERSAL / m_mix
    v1 = R_mix * T1 / p1
    v2 = v1 / eps
    n_p = sum(P_mol.values())
    # pressure after all heat at constant volume (products, v2)
    T_cv = heat_temperature(P_mol, T2, dH, "v")
    p_cv = n_p * R_UNIVERSAL * T_cv / (m_mix * v2)
    if p_cv <= p_max:
        T3 = T4 = T_cv
        p3 = p_cv
        beta = 1.0
    else:
        p3 = p_max
        T3 = p3 * v2 * m_mix / (n_p * R_UNIVERSAL)
        mix = mixture_enthalpy(P_mol, T3) - mixture_enthalpy(P_mol, T2) - n_p * R_UNIVERSAL * (T3 - T2)
        T4 = heat_temperature(P_mol, T3, dH - mix, "p")
        beta = T4 / T3
    alpha = p3 / p2
    w = dual_cycle_work(p1, v1, eps, alpha, beta, n1, n2)
    T5 = T4 * (beta / eps) ** (n2 - 1)
    p5 = p3 * (beta / eps) ** n2
    q = dH / m_mix
    eff = w / q
    return DualCycleResult(rxn, lhv_fuel, q, alpha, beta, {1: T1, 2: T2, 3: T3, 4: T4, 5: T5},
                           {1: p1, 2: p2, 3: p3, 4: p3, 5: p5}, v1, w, w / (v1 - v2), eff,
                           power / (eff * lhv_fuel))
