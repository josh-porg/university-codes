"""Side-by-side results of the AE 571 final-project codes (own and other teams'), heavy gasoline and
hydrogen at ER 0.8, compression ratio 16, peak pressure 5.5 MPa, 55 kW.

The reference row is :func:`unicodes.thermo.engine_cycles.dual_cycle`
(heat of combustion at T2 released into the products, constant volume then
constant pressure, mixture gas constant). The teams differ in where they
evaluate the heating value, whose heat capacities absorb the heat (products
or reactants, mean or local values) and which gas constant sets ``v1``.
"""

import importlib
import sys
from pathlib import Path

from unicodes.thermo.engine_cycles import dual_cycle

sys.path.insert(0, str(Path(__file__).parent))
TEAMS = ["group3", "group7", "group9", "onedrive_team", "team2", "team4", "team5", "team6", "team10"]

for fuel in ("gasoline (heavy)", "hydrogen"):
    print(f"\n{fuel}")
    print(f"{'code':16s} {'LHV kJ/kg':>10s} {'alpha':>7s} {'beta':>7s} {'w kJ/kg':>9s} {'mep kPa':>9s} {'eta':>7s} {'fuel g/s':>9s}")
    ref = dual_cycle(fuel)
    print(f"{'reference':16s} {ref.lhv_mixture / 1e3:10.2f} {ref.alpha:7.4f} {ref.beta:7.4f} {ref.work / 1e3:9.2f} "
          f"{ref.mep / 1e3:9.1f} {ref.efficiency:7.4f} {ref.fuel_flow * 1e3:9.4f}")
    for name in TEAMS:
        mod = importlib.import_module(name)
        run = mod.run_part2 if (name == "team5" and fuel == "hydrogen") else mod.run
        r = run(fuel)
        print(f"{name:16s} {r['lhv'] / 1e3:10.2f} {r['alpha']:7.4f} {r['beta']:7.4f} {r['work'] / 1e3:9.2f} "
              f"{r['mep'] / 1e3:9.1f} {r['eta']:7.4f} {r['fuel_flow'] * 1e3:9.4f}")
