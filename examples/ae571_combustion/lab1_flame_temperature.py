"""AE 571 lab 1: adiabatic flame temperature by Newton iteration (``lab_1_poznanski``,
``getHFromTables``, ``getTable1``, ``getTable2``, ``lab1_hFromTable``, ``getAdiabFlameTemp``,
``str2Elements``, ``string_chem_formula``).

Methane and propane (and hydrogen) burn completely with air or pure oxygen
at equivalence ratios 1.0 to 0.2 from reactants at 300 K. The lab read the
species polynomials from two Excel tables that are not in the repository;
:mod:`unicodes.thermo.combustion` carries the NASA polynomials instead
(diesel, which only had table data, is left out). The lab's lean balance
used ``a = a_stoich * ER``; it is ``a_stoich / ER``.
"""

import matplotlib.pyplot as plt

from unicodes.thermo import adiabatic_flame_temperature, atom_balance, enthalpy_molar, parse_formula
from unicodes.thermo.combustion import cp_molar

print("str2Elements check:", parse_formula("C22H10PuCrN2O5"), parse_formula("CH4O"))
T_amb = 300.0
ratios = [1.0, 0.8, 0.6, 0.4, 0.2]
for fuel in ("methane", "propane", "hydrogen"):
    for oxidizer in ("air", "O2"):
        rows = []
        fig, ax = plt.subplots()
        for er in ratios:
            rxn = atom_balance(fuel, er)
            n2 = rxn.n2 if oxidizer == "air" else 0.0
            H_r = (enthalpy_molar(rxn.fuel.species, T_amb) + rxn.a * enthalpy_molar("O2", T_amb)
                   + n2 * enthalpy_molar("N2", T_amb))

            def H_p(T):
                return (rxn.b * enthalpy_molar("CO2", T) + rxn.c * enthalpy_molar("H2O", T)
                        + rxn.d * enthalpy_molar("O2", T) + n2 * enthalpy_molar("N2", T))

            def dH_p(T):
                return rxn.b * cp_molar("CO2", T) + rxn.c * cp_molar("H2O", T) + rxn.d * cp_molar("O2", T) + n2 * cp_molar("N2", T)

            T, hist = 1000.0, [1000.0]
            for _ in range(50):  # Newton iteration of the lab
                T = T + (H_r - H_p(T)) / dH_p(T)
                hist.append(T)
                if abs((H_r - H_p(T)) / H_r) < 1e-3 * 1e-2:
                    break
            rows.append((er, T, adiabatic_flame_temperature(fuel, er, oxidizer, T_amb), len(hist) - 1))
            ax.plot(hist, "-o", ms=3, label=f"ER {er}")
        ax.set(title=f"{fuel} / {oxidizer}: Newton iterations", xlabel="iteration", ylabel="T (K)")
        ax.legend()
        print(f"\n{fuel} with {oxidizer}:")
        for er, T_newton, T_brent, n in rows:
            print(f"  ER {er:.1f}: T_ad = {T_newton:7.1f} K ({n} Newton steps; brentq {T_brent:7.1f} K)")
plt.show()
