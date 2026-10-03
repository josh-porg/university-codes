"""Combustion of hydrocarbon fuels with air: atom balance and species thermo.

Ported from the AE 571 final project (team 8). Molar units throughout:
enthalpy in J/mol (= kJ/kmol), cp in J/(mol K).

Corrections relative to the MATLAB:

* NASA polynomial enthalpy used ``a5/5*T`` instead of ``a5/5*T**5``.
* Fuel cp used ``A5/theta`` instead of ``A5/theta**2``, and the heavy
  gasoline enthalpy dropped the ``A6`` term (Heywood, Table 4.8).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

R_UNIVERSAL = 8.314  # J/(mol K)
AIR_N2_PER_O2 = 3.76

# Old-format NASA 7-coefficient polynomials (first six coefficients):
# (high: T >= 1000 K, low: T < 1000 K)
NASA7 = {
    "CO2": (
        (0.04453623e2, 0.03140168e-1, -0.12784105e-5, 0.02393996e-8, -0.16690333e-13, -0.04896696e6),
        (0.02275724e2, 0.09922072e-1, -0.10409113e-4, 0.06866686e-7, -0.0211728e-10, -0.04837314e6),
    ),
    "H2O": (
        (0.02672145e2, 0.03056293e-1, -0.0873026e-5, 0.12009964e-9, -0.06391618e-13, -0.02989921e6),
        (0.03386842e2, 0.03474982e-1, -0.06354696e-4, 0.06968581e-7, -0.02506588e-10, -0.03020811e6),
    ),
    "O2": (
        (0.03697578e2, 0.06135197e-2, -0.12588420e-6, 0.01775281e-9, -0.11364354e-14, -0.12339301e4),
        (0.03212936e2, 0.11274864e-2, -0.0575615e-5, 0.13138773e-8, -0.08768554e-11, -0.1005249e4),
    ),
    "N2": (
        (0.0292664e2, 0.14879768e-2, -0.0568476e-5, 0.10097038e-9, -0.06753351e-13, -0.09227977e4),
        (0.03298677e2, 0.14082404e-2, -0.03963222e-4, 0.05641515e-7, -0.02444854e-10, -0.10208999e4),
    ),
    "H2": (
        (0.02991423e2, 0.07000644e-2, -0.05633828e-6, -0.09231578e-10, 0.15827519e-14, -0.0835034e4),
        (0.03298124e2, 0.08249441e-2, -0.08143015e-5, -0.09475434e-9, 0.04134872e-11, -0.10125209e4),
    ),
}

# Heywood, Internal Combustion Engine Fundamentals, Table 4.8:
# cp = 4.184 (A1 + A2 t + A3 t^2 + A4 t^3 + A5 / t^2),  t = T / 1000
# h  = 4184 (A1 t + A2 t^2/2 + A3 t^3/3 + A4 t^4/4 - A5 / t + A6)
HEYWOOD_FUELS = {
    "C8.26H15.5": (-24.078, 256.63, -201.68, 64.75, 0.5808, -27.562),
    "C7.76H13.1": (-22.501, 227.99, -177.26, 56.048, 0.4845, -17.578),
}


@dataclass(frozen=True)
class Fuel:
    """A CxHy fuel."""

    name: str
    x: float  # carbon atoms
    y: float  # hydrogen atoms
    molar_mass: float  # g/mol
    species: str | None = None  # key for thermo data, if available


FUELS = {
    "gasoline (heavy)": Fuel("Gasoline (Heavy)", 8.26, 15.5, 114.8, "C8.26H15.5"),
    "gasoline (light)": Fuel("Gasoline (Light)", 7.76, 13.1, 106.4, "C7.76H13.1"),
    "hydrogen": Fuel("Hydrogen", 0.0, 2.0, 2.0, "H2"),
}


@dataclass(frozen=True)
class LeanCombustion:
    """Coefficients of  CxHy + a (O2 + 3.76 N2) -> b CO2 + c H2O + d O2 + 3.76 a N2."""

    fuel: Fuel
    equivalence_ratio: float
    a: float
    b: float
    c: float
    d: float

    @property
    def n2(self) -> float:
        return AIR_N2_PER_O2 * self.a


def atom_balance(fuel: str | Fuel, equivalence_ratio: float) -> LeanCombustion:
    """Balance complete lean (or stoichiometric) combustion of one mole of fuel in air."""
    if isinstance(fuel, str):
        fuel = FUELS[fuel.lower()]
    if not 0 < equivalence_ratio <= 1:
        raise ValueError("Only lean or stoichiometric mixtures (0 < ER <= 1) are supported")
    a = (fuel.x + fuel.y / 4) / equivalence_ratio
    b = fuel.x
    c = fuel.y / 2
    d = a - b - c / 2
    return LeanCombustion(fuel, equivalence_ratio, a, b, c, d)


def _nasa_coeffs(species: str, T: float):
    high, low = NASA7[species]
    return np.array(high if T >= 1000 else low)


def cp_molar(species: str, T: float) -> float:
    """Constant-pressure molar specific heat, J/(mol K)."""
    if species in HEYWOOD_FUELS:
        A1, A2, A3, A4, A5, _ = HEYWOOD_FUELS[species]
        t = T / 1000
        return 4.184 * (A1 + A2 * t + A3 * t**2 + A4 * t**3 + A5 / t**2)
    a = _nasa_coeffs(species, T)
    return R_UNIVERSAL * (a[0] + a[1] * T + a[2] * T**2 + a[3] * T**3 + a[4] * T**4)


def enthalpy_molar(species: str, T: float) -> float:
    """Absolute (formation + sensible) molar enthalpy, J/mol."""
    if species in HEYWOOD_FUELS:
        A1, A2, A3, A4, A5, A6 = HEYWOOD_FUELS[species]
        t = T / 1000
        return 4184 * (A1 * t + A2 * t**2 / 2 + A3 * t**3 / 3 + A4 * t**4 / 4 - A5 / t + A6)
    a = _nasa_coeffs(species, T)
    return R_UNIVERSAL * (
        a[0] * T + a[1] / 2 * T**2 + a[2] / 3 * T**3 + a[3] / 4 * T**4 + a[4] / 5 * T**5 + a[5]
    )


def heat_of_combustion(rxn: LeanCombustion, T: float, T_ref: float = 298.0) -> float:
    """Lower heating value released per mole of fuel at temperature ``T``, J/mol.

    Reactant and product enthalpies are taken at ``T_ref`` and carried to
    ``T`` with cp averaged between the two temperatures, as in the project.
    """
    sp = rxn.fuel.species

    def cp_avg(s):
        return (cp_molar(s, T) + cp_molar(s, T_ref)) / 2

    dT = T - T_ref
    H_r = (enthalpy_molar(sp, T_ref) + rxn.a * enthalpy_molar("O2", T_ref) + rxn.n2 * enthalpy_molar("N2", T_ref)) + (
        cp_avg(sp) + rxn.a * cp_avg("O2") + rxn.n2 * cp_avg("N2")
    ) * dT
    H_p = (
        rxn.b * enthalpy_molar("CO2", T_ref)
        + rxn.c * enthalpy_molar("H2O", T_ref)
        + rxn.d * enthalpy_molar("O2", T_ref)
        + rxn.n2 * enthalpy_molar("N2", T_ref)
    ) + (rxn.b * cp_avg("CO2") + rxn.c * cp_avg("H2O") + rxn.d * cp_avg("O2") + rxn.n2 * cp_avg("N2")) * dT
    return H_r - H_p
