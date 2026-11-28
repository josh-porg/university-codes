"""AE 573 homework 10 (``HW_10``): mixed-flow afterburning turbofan.

M0 = 2, p0 = 10 kPa, T0 = 223 K; pi_d = 1 (the file first set 0.9, then
overwrote it), pi_f = 1.9 (e_f = 0.9), e_c = 0.9, bypass duct pi_fd = 0.99,
Tt4 = 1600 K, pi_b = 0.95 (eta_b = 0.98, Qr = 42 MJ/kg), e_t = 0.8,
eta_m = 1 (first typed as 95, then overwritten), core Mach number at the
mixer M5 = 0.5, mixer friction loss 0.98; afterburner Tt7 = 2000 K,
pi_AB = 0.92, eta_AB = 0.98, gamma_AB = 1.3, cp_AB = 1241; pi_n = 0.95 and
the nozzle exit pressure p9 = 3.8 p0 (``pi_inlet_outlet = p_9/p_0``).
Turbine gases gamma 1.33, cp 1152.

The file never gave the overall compressor pressure ratio, which the
relations need: with it missing, the relation solver could only find
quantities upstream of the compressor. It is a command-line argument here
(default 15, as in the final exam engine with the same fan). The bypass
ratio follows from matching the mixer total pressures. A pi_c sweep
shows how alpha, the specific thrust and the fuel consumption move::

    python hw_10.py --pi-c 20
"""

import argparse

from unicodes.propulsion_cycles import AIR, Gas, mixed_turbofan

p = argparse.ArgumentParser()
p.add_argument("--pi-c", type=float, default=15.0)
p.add_argument("--p9-p0", type=float, default=3.8, help="nozzle exit over ambient pressure")
a = p.parse_args()

hot, ab = Gas(1.33, 1152.0), Gas(1.3, 1241.0)
kw = dict(pi_f=1.9, Tt4=1600, Qr=42e6, M6=0.5, gas_c=AIR, gas_t=hot, pi_d=1.0, pi_b=0.95, pi_fd=0.99,
          pi_m_friction=0.98, pi_n=0.95, e_c=0.9, e_f=0.9, e_t=0.8, eta_b=0.98, eta_m=1.0, Tt7=2000, gas_ab=ab,
          pi_ab=0.92, eta_ab=0.98, p9=a.p9_p0 * 10e3)
r = mixed_turbofan(2.0, 223.0, 10e3, pi_c=a.pi_c, **kw)
print(f"pi_c = {a.pi_c}")
print(r.summary())

print("\n pi_c   alpha   F/mdot0 [N s/kg]   TSFC [mg/(N s)]")
for pi_c in (8, 10, 12, 15, 20, 25, 30):
    try:
        s = mixed_turbofan(2.0, 223.0, 10e3, pi_c=pi_c, **kw)
        print(f"{pi_c:5.0f} {s.alpha:7.3f} {s.specific_thrust:14.2f} {s.tsfc * 1e6:17.3f}")
    except ValueError as err:
        print(f"{pi_c:5.0f}  {err}")
