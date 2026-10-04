"""AE 571 final project: limited-pressure (dual) cycle of a spark-ignition engine burning heavy
gasoline (``poznanski_main``, ``team_code``, ``AE571_FinalProject_Team8``).

Compression ratio 16, ER 0.8, polytropic compression (n = 1.38) and expansion
(n = 1.25), ambient 100 kPa / 298 K, peak pressure 5.5 MPa (the team code
asked for it with ``input``; ``--p-max``). The heat of combustion comes from
:func:`unicodes.thermo.heat_of_combustion` at the end-of-compression
temperature, as in the team code.

Fixes: the team code computed ``T4 = (T3 + q - cv) / cp`` (mixing a
temperature and energies) and ``cv = cp * R``; the constant-pressure
burning is closed here with ``T4 = T3 + q_remaining / cp`` and
``cv = cp - R``. ``poznanski_main`` used the molar-mass sum of products in
place of the mixture molar mass for mass fractions; mole-weighted values are used.
"""

import argparse

import numpy as np

from unicodes.thermo import atom_balance, heat_of_combustion
from unicodes.thermo.combustion import cp_molar

p = argparse.ArgumentParser()
p.add_argument("--p-max", type=float, default=5.5e6)
p.add_argument("--er", type=float, default=0.8)
p.add_argument("--eps", type=float, default=16.0)
a = p.parse_args()

n1, n2_exp = 1.38, 1.25
p1, T1 = 100e3, 298.0
rxn = atom_balance("gasoline (heavy)", a.er)
p2 = p1 * a.eps**n1
T2 = T1 * a.eps ** (n1 - 1)
dHc = heat_of_combustion(rxn, T2)  # J per mol fuel
M_fuel = rxn.fuel.molar_mass / 1000
m_mix = M_fuel + rxn.a * 0.032 + rxn.n2 * 0.028  # kg per mol fuel
q = dHc / m_mix  # J per kg mixture
moles = np.array([rxn.b, rxn.c, rxn.d, rxn.n2])
M = np.array([44.01, 18.02, 32.0, 28.01])
x = moles / moles.sum()
print(f"LHV {dHc / M_fuel / 1e6:.2f} MJ/kg fuel, heat release {q / 1e6:.3f} MJ/kg mixture")
print("product mole fractions [CO2, H2O, O2, N2]:", np.round(x, 4), " mass fractions:", np.round(x * M / (x @ M), 4))

cp_mix = (cp_molar(rxn.fuel.species, T2) + rxn.a * cp_molar("O2", T2) + rxn.n2 * cp_molar("N2", T2)) / m_mix
R = 287.05
cv = cp_mix - R
v1 = R * T1 / p1
v2 = v1 / a.eps
alpha = a.p_max / p2  # constant-volume pressure rise
T3 = alpha * T2
q_v = cv * (T3 - T2)
if q_v > q:
    raise SystemExit("peak pressure too high for the available heat release")
T4 = T3 + (q - q_v) / cp_mix  # constant-pressure burning
beta = T4 / T3  # cut-off ratio
delta = a.eps / beta  # expansion ratio
W = p1 * v1 * a.eps ** (n1 - 1) * (alpha * (beta - 1) + alpha * beta / (n2_exp - 1) * (1 - delta ** (1 - n2_exp))
                                  - (1 - a.eps ** (1 - n1)) / (n1 - 1))
print(f"T2 = {T2:.1f} K, p2 = {p2 / 1e6:.3f} MPa, alpha = {alpha:.3f}, T3 = {T3:.1f} K, T4 = {T4:.1f} K, beta = {beta:.3f}")
print(f"net work {W / 1e3:.1f} kJ/kg, MEP {W / (v1 - v2) / 1e6:.3f} MPa, thermal efficiency {W / q:.3f}")
