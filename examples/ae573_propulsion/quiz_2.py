"""AE 573 quiz 2 (``Quiz_2``): afterburning turbojet.

M0 = 2, p0 = 25 kPa, T0 = 228.15 K, mdot = 100 kg/s; pi_d = 0.9, pi_c = 15
(e_c = 0.9), pi_b = 0.95 (eta_b = 0.99, Qr = 42 MJ/kg), Tt4 = 1760 K,
e_t = 0.85, pi_AB = 0.98, pi_n = 0.95; the nozzle exit is M9 = 2.62,
T9 = 603 K. Cold gas gamma 1.4, hot gas gamma 1.33, cp from gamma and
R = 287.

Corrections to the MATLAB relations: ``tau_c_lp = Tt_2 / Tt_2_5`` and
``tau_5 = Tt_4 / Tt_5`` were inverted, ``pi_r = pt_0 / P_0`` used an
undefined ``P_0``, and the burner and afterburner heat balances shared one
``Qr`` relation with the wrong mass flows. The file also gave
``pi_t = 0.936``, which cannot hold together with the rest: with
``tau_t = pi_t^((gamma_t-1) e_t / gamma_t) = 0.987`` the turbine could only
drive the compressor with Tt4 above 30 000 K. So Tt4 = 1760 K is taken as
given and the turbine follows from the shaft balance (eta_m = 1) and its
polytropic efficiency. The given nozzle exit fixes Tt7, the afterburner fuel
and the exit pressure (the nozzle is not fully expanded, so the thrust
carries a pressure term when it is not fully expanded). The quiz exit state
comes out at Tt9 = Tt5 and p9 = p0, i.e. it describes the engine with the
afterburner off, which checks the turbine and burner values above; the
afterburner fuel is then zero. The file's ``Fn = 4588.86`` is printed for
comparison.
"""

import numpy as np

from unicodes.propulsion_cycles import Gas, burner_fuel_air_ratio, compressor, freestream

R = 287.0
c, t = Gas(1.4, 1.4 * R / 0.4), Gas(1.33, 1.33 * R / 0.33)
M0, p0, T0, mdot = 2.0, 25e3, 228.15, 100.0
pi_d, pi_c, e_c, pi_b, eta_b, Qr, Tt4 = 0.9, 15.0, 0.9, 0.95, 0.99, 42e6, 1760.0
e_t, pi_AB, pi_n = 0.85, 0.98, 0.95
M9, T9 = 2.62, 603.0

V0, st0, tau_r, pi_r = freestream(M0, T0, p0, c)
pt2, Tt2 = st0.pt * pi_d, st0.Tt
tau_c, eta_c = compressor(pi_c, e_c, c)
Tt3, pt3 = Tt2 * tau_c, pt2 * pi_c
f = burner_fuel_air_ratio(Tt3, Tt4, Qr, eta_b, c, t)
Tt5 = Tt4 - c.cp * (Tt3 - Tt2) / ((1 + f) * t.cp)
tau_t = Tt5 / Tt4
pi_t = tau_t ** (t.gamma / ((t.gamma - 1) * e_t))
pt5 = pt3 * pi_b * pi_t
pt9 = pt5 * pi_AB * pi_n
ratio9 = 1 + (t.gamma - 1) / 2 * M9**2
Tt9 = T9 * ratio9
p9 = pt9 / ratio9 ** (t.gamma / (t.gamma - 1))
f_AB = (1 + f) * t.cp * (Tt9 - Tt5) / (eta_b * Qr - t.cp * Tt9)
V9 = M9 * np.sqrt(t.gamma * t.R * T9)
m9 = 1 + f + f_AB
F = m9 * V9 - V0 + m9 * t.R * T9 / V9 * (1 - p0 / p9)

print(f"V0 = {V0:.2f} m/s, Tt2 = {Tt2:.2f} K, pt2 = {pt2:.0f} Pa")
print(f"Tt3 = {Tt3:.2f} K (eta_c = {eta_c:.4f}), pt3 = {pt3:.0f} Pa")
print(f"Tt4 = {Tt4:.1f} K, f = {f:.5f}; turbine tau_t = {tau_t:.4f}, pi_t = {pi_t:.4f}; Tt5 = {Tt5:.1f} K, "
      f"pt5 = {pt5:.0f} Pa")
print(f"nozzle: pt9 = {pt9:.0f} Pa, p9 = {p9:.0f} Pa (p0 = {p0:.0f}), Tt9 = Tt7 = {Tt9:.1f} K, V9 = {V9:.1f} m/s")
print(f"f_AB = {f_AB:.5f}, total fuel {mdot * (f + f_AB):.3f} kg/s")
print(f"F/mdot0 = {F:.2f} N s/kg (file: Fn = 4588.86), F = {F * mdot / 1e3:.2f} kN, "
      f"TSFC = {(f + f_AB) / F * 1e6:.3f} mg/(N s)")
